#' A new quasi-likelihood method to analyze the RNA-Seq data
#' @param data.norm Raw count data after normalization
#' @param design1 Design matrix for the full model
#' @param design2 Design matrix for the redcued model
#' @param splmodel The average variance function estimated from the data (spline)
#' @return Returns raw p-values
#' @examples
#' ## Load example data
#'
#' library(airway)
#' library(edgeR)
#' data(airway)
#' airraw <- assays(airway)$count
#' air <- airraw[,c(1,3,5,7,2,4,6,8)]
#' n1 <- 4
#' n2 <- 4
#' n <- 8
#' trt<-c(rep(1,n1),rep(2,n2))
#' ## Filtering
#' IQR <- apply(air,1,IQR)
#' airf <- air[IQR != 0 & rowMeans(air)>1,]
#'
#' ## Normalization
#'
#' libsize = apply(airf[,1:n],2,sum)
#' nrmfactor <-calcNormFactors(airf[,1:n],method = "TMM")*libsize
#' air.norm <- as.data.frame(t(t(airf[,1:n])/nrmfactor))*mean(nrmfactor)
#'
#' ## Estimate average variance function
#'
#' air.norm$mean <- apply(air.norm[,1:n],1,mean)
#' air.norm$mean1 <- apply(air.norm[,1:n1],1,mean)
#' air.norm$mean2 <- apply(air.norm[,(n1+1):n],1,mean)
#' air.norm$var <- apply(air.norm[,1:n],1,var)
#' air.norm$var1 <- apply(air.norm[,1:n1],1,var)
#' air.norm$var2 <- apply(air.norm[,(n1+1):n],1,var)
#' meanest <- c(air.norm$mean1,air.norm$mean2)
#' varest <- c(air.norm$var1,air.norm$var2)
#' modnorm <- smooth.spline(meanest[varest != 0],log(varest[varest != 0]))
#'
#' ## Analyze the normalized data
#' 
#' ## Use a subset for testing
#' air.norm <- air.norm[1:3000,]
#' dsgn <- model.matrix(~as.factor(trt))
#' air.norm$QuasiDE.rawp <- QuasiDE(air.norm[,1:n],design1=dsgn,splmodel=modnorm)
#' air.norm$QuasiDE.adjp <- p.adjust(air.norm$QuasiDE.rawp,method="BH")
#' sum(air.norm$QuasiDE.adjp<0.05)
#' @export
QuasiDE <- function(data.norm,design1,design2=NULL,splmodel){
  j <- 1
  res <- matrix(NA,nrow=nrow(data.norm),ncol=2)
  repeat{
    yobs=as.numeric(data.norm[j,])
    stf1 <- as.data.frame(cbind(yobs,design1))
    sfit1 <- glm(yobs~.,data=stf1)
    stf1$u <- predict(sfit1)
    if (! is.null(design2)){
      stf2 <-  as.data.frame(cbind(yobs,design2))
      sfit2 <- glm(yobs~.,data=stf2)
      stf2$u <- predict(sfit2)
    }
    fit1 <- tryCatch({glm(yobs~design1, family = quasi.spl(variance = "spline", link = "log",splmodel=splmodel))}, error = function(e){} )
    if (is.null(fit1)) fit1 <- tryCatch({glm(yobs~design1, family = quasi.spl(variance = "spline", link = "log",splmodel=splmodel),mustart=stf1$u)}, error = function(e){} )
    if (is.null(fit1)) fit1 <- sfit1
    if (! is.null(design2)){
      fit2 <- tryCatch({glm(yobs~design2, family = quasi.spl(variance = "spline", link = "log",splmodel=splmodel))}, error = function(e){} )
      if (is.null(fit2)) fit2 <- tryCatch({glm(yobs~design2, family = quasi.spl(variance = "spline", link = "log",splmodel=splmodel),mustart=stf2$u)}, error = function(e){} )
      if (is.null(fit2)) fit2 <- sfit2
    }

    if ((!is.null(fit1)) & (!is.null(design2))) {
      if (!is.null(fit2)) {
        res[j,1] <- (fit2$deviance-fit1$deviance)/(ncol(design1)-ncol(design2))
        res[j,2] <- sum(residuals(fit1,"pearson")^2)/fit1$df.residual
        num.df <- (ncol(design1)-ncol(design2))
        den.df <- fit1$df.residual
      }
    }
    if (is.null(design2) & (!is.null(fit1))) {
      res[j,1] <- (fit1$null.deviance-fit1$deviance)/(ncol(design1)-1)
      res[j,2] <- sum(residuals(fit1,"pearson")^2)/fit1$df.residual
      num.df <- (ncol(design1)-1)
      den.df <- fit1$df.residual
    }
    j <- j + 1
    if ((j %% 1000)==0) print(paste(j,"genes evaluated ..."))
    if (j>nrow(data.norm)) break
  }
  LRT <- res[,1]
  phi <- res[,2]
  logphi <- log(phi)
  logphi[logphi == -Inf] <- min(logphi[logphi != -Inf])
  logphi[logphi == Inf] <- max(logphi[logphi != Inf])
  naind <- is.na(logphi)
  logphi <- logphi[! naind]
  data.mn <- apply(data.norm,1,mean)[! naind]
  spline.fit <- smooth.spline(log(data.mn), logphi)
  y2 <- phi/exp(predict(spline.fit,log(data.mn))$y)
  phi.hat <- y2
  phi.hat[phi.hat <= 0] <- min(phi.hat[phi.hat > 0])
  z <- log(phi.hat)
  z[z == Inf] <- max(z[z != Inf])
  z[z == -Inf] <- min(z[z != -Inf])
  mnz <- mean(z)
  d0arg <- var(z) - trigamma((den.df)/2)
  if (d0arg > 0) {
    dif <- function(x, y) abs(trigamma(x) - y)
    inverse.trigamma <- function(y) optimize(dif, interval = c(0,10000), y = y)$minimum
    d0 <- 2 * inverse.trigamma(d0arg)
    phi0 <- exp(mnz - digamma((den.df)/2) + digamma(d0/2) - log(d0/(den.df)))
    phi.shrink <- ((den.df) * phi.hat + d0 * phi0)/(den.df + d0)
  }
  else {
    phi.shrink <- rep(exp(mnz), length(z))
    d0 <- Inf
    phi0 <- exp(mnz)
  }
  if (d0== Inf) {
    phi.spline <- exp(predict(spline.fit,log(data.mn))$y)
  } else {
    phi.spline <- (d0 * exp(predict(spline.fit,log(data.mn))$y) * phi0 + den.df * phi)/(d0 + den.df)
  }
  rawp <- 1-pf(LRT/phi.spline, df1 = num.df, df2=d0 + den.df)
  return(rawp)
}

#' Modified Quasi-likelihood function used in GLM function
#' @param link link function
#' @param variance variance function name
#' @param splmodel the spline object used as variance function
#' @return Returns an object for the QuasiDE function to fit GLM model using qusi-likelihood method

quasi.spl <- function (link = "log", variance = "spline", splmodel)
{
  # trapz - trapezoid rule to calculate Area under the curve
  trapz <- function (x,y) {
    area <- NULL
    if (length(x) == 1) area <- 0
    else area <- sum((y[-1]+y[-length(y)])/2*(x[-1]-x[-length(x)]))
    return(area)
  }
  areadev.unit <- function(y,mu,wt){
    res <- NULL
    if ((length(mu) == 1) && length(mu)<length(y)) m = rep(mu,length(y))
    else m = mu
    for (i in 1:length(y)){
      x <- seq(m[i], y[i],len = 501)
      ql <- 2*(y[i]-x)/(exp(predict(splmodel,x)$y))
      if (sum(is.nan(ql))>0) {
        x <- x[! is.nan(ql)]
        ql <- ql[! is.nan(ql)]
      }
      res <- c(res,trapz(x,ql))
    }
    return(res)
  }
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "identity",
                      "inverse", "log", "1/mu^2", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    stats <- link
    linktemp <- if (!is.null(stats$name))
      stats$name
    else deparse(linktemp)
  }
  vtemp <- substitute(variance)
  if (!is.character(vtemp))
    vtemp <- deparse(vtemp)
  variance_nm <- vtemp
  switch(vtemp, `spline` = {
    varfun <- function(mu) exp(predict(splmodel,mu)$y)
    validmu <- function(mu) all(mu > 0)
    dev.resids <- function(y, mu,wt) areadev.unit(y,mu,wt)
    initialize <- expression({
      n <- rep.int(1, nobs)
      mustart <- y + 0.1 * (y == 0)
    })
  }, variance_nm <- NA)
  if (is.na(variance_nm)) {
    if (is.character(variance))
      stop(gettextf("'variance' \"%s\" is invalid: possible values are \"spline\"",
                    variance_nm), domain = NA)
    varfun <- variance$varfun
    validmu <- variance$validmu
    dev.resids <- variance$dev.resids
    initialize <- variance$initialize
    variance_nm <- variance$name
  }
  aic <- function(y, n, mu, wt, dev) NA
  structure(list(family = "quasi", link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = varfun, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, varfun = variance_nm),
            class = "family")
}
