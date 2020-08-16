
.get.S <- function(t, time, surv)	# assumes that t is a scalar
{
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(1)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(surv[n])
  #if (t>=time[1] && t<=time[n])
  return( surv[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}


#' Diagnostics of non-overlap in the TVW estimator for a binary covariate
#'
#'
#' Visual test for the marginal (i.e. for a single covariate) violation of
#' positivity (overlap).
#' The diagnostics are done by visual comparing between the proportions of 1
#' in a binary covariate \code{Z} corresponding to
#' observations at risk at each event time, for TVW versus CW estimators.
#' This test can identify practical, TVW-specific violation of positivity.
#' But it cannot test theoretical violation of positivity.
#'
#'
#' @param Z a binary covariate for which we want to test whether the TVW-specific positivity is not violated.
#' @param X a vector of observed (possibly right-censored) lifetimes.
#' @param T a vector of truncation times.
#' @param delta a censoring indicator. It is assumed that 0 corresponds to censoring.
#' @param censoring censoring can get one of the three values "after" (default), "before", "none".
#'
#' @return The function plots the time-varying proportions of 1 and does not return a value.
#'
#' @examples
#' # under construction
#'
#' @export
positivity.binary <- function(Z, X, T, delta, ylim=NA, xlim=NA)
{
  nev <- sum(delta)
  X.m <- matrix(X[delta==1], n, nev, byrow=TRUE)
  X.i <- matrix(X, n, nev, byrow=FALSE)
  Z1 <- matrix(Z, n, nev, byrow=FALSE)

  d.Z.cw <- data.frame(X=X.m[X.i>=X.m], Z=Z1[X.i>=X.m], type=1)
  cw.avg = aggregate(list(Z = d.Z.cw$Z), list(OS = factor(d.Z.cw$X)), mean)


  # TVW:
  T.i <- matrix(T, n, nev, byrow=FALSE)
  d.Z.tvw <- data.frame(X=X.m[(X.i>=X.m)  & (X.m >=T.i)],
                        Z=Z1[(X.i>=X.m)  & (X.m >=T.i)], type=2)
  tvw.avg = aggregate(list(Z = d.Z.tvw$Z), list(OS = factor(d.Z.tvw$X)), mean)

  d.Z <- rbind(cbind(cw.avg, type=1), cbind(tvw.avg, type=2))
  #d.Z <- d.Z[d.Z$X %in% X[ch],]
  d.Z$OS <- as.numeric(paste(d.Z$OS))
  d.Z$OS <- round(d.Z$OS, dig=2)
  #d.Z$OS <- as.factor(d.Z$OS)
  d.Z$type <- factor(d.Z$type)
  levels(d.Z$type) <- c("CW", "TVW")

  # distribution of Z_j in risk sets of CW and TVW estiamtors, compared to Z_i
  p <- ggplot(d.Z, aes(x=OS, y=Z, group=type, col=type, fill=type)) +
    geom_line(aes(linetype=type), size=1.2)+ geom_point(aes(shape=type), size=2.5)+
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold")) + ylim(0,1)+
    xlab(expression(paste(X[i],'*', sep="")))+
    ylab(paste("proportion of ones among survivors to i^th time", sep=""))+
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")+
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold"))
  if (is.na(ylim) && is.na(xlim))  plot(p)
  if (!is.na(ylim))  plot(p + ylim(ylim[1], ylim[2]))
  if (!is.na(xlim))  plot(p + xlim(xlim[1], xlim[2]))
  if (!is.na(ylim) && !is.na(xlim))  plot(p + ylim(ylim[1], ylim[2]) + xlim(xlim[1], xlim[2]))
}

#' Diagnostics of non-overlap in the TVW estimator for a continuous covariate
#'
#'
#' Visual test for the marginal (i.e. for a single covariate) violation of
#' positivity (overlap).
#' The diagnostics are done by visual comparing between the distributions
#' of a continuous covariate corresponding to
#' observations at risk at each event time, for TVW versus CW estimators.
#' This test can identify practical, TVW-specific violation of positivity.
#' But it can neither test theoretical violation of positivity, nor
#' the non-positivity with respect to the lifetime \code{X}  itself.
#'
#'
#' @param Z a continuous covariate for which we want to test whether the TVW-specific positivity is not violated.
#' @param X a vector of observed (possibly right-censored) lifetimes.
#' @param T a vector of truncation times.
#' @param delta a censoring indicator. It is assumed that 0 corresponds to censoring.
#' @param censoring censoring can get one of the three values "after" (default), "before", "none".
#'
#' @return The function plots the distributions and doe not return a value.
#'
#' @examples
#' # under construction
#'
#' @export

positivity.cont <- function(Z, X, T, delta, chosen=NA, ylim=NA, xlim=NA)
{
  nev <- sum(delta)
  if (is.na(chosen[1]))
  {

    if (nev <100) ch <- 1:nev
    else
    {
      by <- round(nev/100, dig=0)
      ch <- seq(1,nev, by=by) # a little more than 100...
    }
  }
  else
  {
    ch <- chosen[1:min(c(100, length(chosen)))]
  }

  X.m <- matrix(X[delta==1], n, nev, byrow=TRUE)
  X.i <- matrix(X, n, nev, byrow=FALSE)
  Z1 <- matrix(Z, n, nev, byrow=FALSE)

  d.Z.cw <- data.frame(X=X.m[X.i>=X.m], Z=Z1[X.i>=X.m], type=1)

  # TVW:
  T.i <- matrix(T.s, n, nev, byrow=FALSE)
  d.Z.tvw <- data.frame(X=X.m[(X.i>=X.m)  & (X.m >=T.i)],
                        Z=Z1[(X.i>=X.m)  & (X.m >=T.i)], type=2)
  d.Z <- rbind(d.Z.cw, d.Z.tvw)
  d.Z <- d.Z[d.Z$X %in% X[ch] ,]
  d.Z$X <- round(d.Z$X, dig=2)
  d.Z$type <- factor(d.Z$type)
  levels(d.Z$type) <- c("CW", "TVW")

  # distribution of Z_j in risk sets of CW and TVW estiamtors, compared to Z_i
  p <- ggplot(d.Z, aes(x=factor(X), y=Z, col=type, fill=type)) +
    geom_boxplot(width=1, alpha=0.6,  position=position_dodge(width=1)) +
    stat_summary(fun=mean, color="black", geom="point",
                 shape=16, size=2, show.legend = FALSE, position=position_dodge(width=1), alpha=0.5)+
    theme(axis.text.x = element_text(face="bold", angle=90),
          axis.text.y = element_text(face="bold")) +# ylim(0,12)+
    xlab(expression(paste(X[i],'*', sep="")))+
    ylab("distribution of Z at risk")+
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")+
    theme(axis.text.x = element_text(face="bold", angle=90),
          axis.text.y = element_text(face="bold"))
  if (is.na(ylim) && is.na(xlim))  plot(p)
  if (!is.na(ylim))  plot(p + ylim(ylim[1], ylim[2]))
  if (!is.na(xlim))  plot(p + xlim(xlim[1], xlim[2]))
  if (!is.na(ylim) && !is.na(xlim))  plot(p + ylim(ylim[1], ylim[2]) + xlim(xlim[1], xlim[2]))
}


# inner function
# if censoring can happen only after sampling:
.get.CW.res.cens <- function(X, T, delta, Z)
{
  f.C <- survfit(Surv(X-T, 1-delta) ~ 1)
  S_c <- sapply(X[delta==1]-T[delta==1], .get.S, time=summary(f.C)$time, surv=summary(f.C)$surv, simplify=TRUE, USE.NAMES=FALSE)

  # using only uncensored observations estimate truncation distribution given Z
  nn <- sum(delta)
  del= rep(1, nn)
  m <- max(c(X[delta==1], T[delta==1]))
  srv.trunc <- Surv(time=m-X[delta==1], time2=m-T[delta==1], event=del)
  #   sol <- try(coxph(srv.trunc ~ Z[delta==1,], weights=1/S_c) ) # if Z is a matrix....
  sol <- try(coxph(srv.trunc ~ Z[delta==1], weights=1/S_c) )
  if (class(sol)!="coxph")
  {
    cat("coxph: Error occurred in 'get.CW.res.cens'.","\n")
  }
  alpha <- matrix(sol$coef, length(sol$coef), 1)

  bh <- basehaz(sol, centered=FALSE)
  H <- bh$hazard # ordered cumulative hazard
  K.S.NA <- exp(-H) # in the order of m-T[n], m-T[n-1], ...
  n1 <- length(bh$hazard)

  K.CDF.forw <- K.S.NA[n1:1] # this is only baseline pr-ties for z=0!!!!
  base <- 1/(1-sapply(X[delta==1], .get.S, time=m-bh$time[n1:1], surv=1-K.CDF.forw, simplify=TRUE, USE.NAMES=FALSE))

  #  s <- as.vector(exp(Z[delta==1,] %*% alpha) )
  s <- as.vector(exp(Z[delta==1] * c(alpha)) )
  w <- base^s/S_c

  f <- survfit(Surv(X[delta==1], delta[delta==1]) ~ 1, weights=w)
  list(time=summary(f)$time, surv=summary(f)$surv)
}

# inner function for CW, allows ties:
# if censoring can happen before sampling:
.get.CW <- function(X, T, delta, Z)
{
  n <- length(T)
  m <- max(c(X, T))
  d.t <- rep(1,n)
  sol <- try(coxph(Surv(time=m-X, time2=m-T, event=d.t) ~ Z) )
  if (class(sol)!="coxph")
  {
    cat("coxph: Error occurred in coxph(srv.trunc ~ Z).","\n")
  }
  alpha <- matrix(sol$coef, length(sol$coef), 1)
  bh <- basehaz(sol, centered=FALSE)
  H <- bh$hazard # ordered cumulative hazard
  K.S.NA <- exp(-H) # in the order of m-T[n], m-T[n-1], ...
  n1 <- length(bh$hazard)

  K.CDF.forw <- K.S.NA[n1:1] # this is only baseline pr-ties for z=0!!!!
  base <- 1/(1-sapply(X, .get.S, time=m-bh$time[n1:1], surv=1-K.CDF.forw, simplify=TRUE, USE.NAMES=FALSE))

  s <- as.vector(exp(Z %*% alpha ) )
  w <- base^s #

  f <- survfit(Surv(X, delta) ~ 1, weights=w) # survfit treats the ties in X correctly, but smth was wrong there in the length times
  list(time=summary(f)$time, surv=summary(f)$surv)
}



# assumes Z is a vector of one covariate !!!!!!!!!!!!!!

#' Case-weight estimation of survival probability
#'
#'
#' The estimator assumes the Cox proportional hazards model for truncation time given covariates, which is
#' used to obtain "case weights" and reweigh the contribution of the sampled observations accordingly.
#' For further details we refer to the paper of Vakulenko-Lagun et al. (2020).
#'
#'
#' @param X a vector of observed (possibly right-censored) lifetimes.
#' @param T a vector of truncation times.
#' @param delta a censoring indicator. It is assumed that 0 corresponds to censoring.
#' @param Z a matrix of covariates (common causes of \code{T} and \code{X}) that will be used for estimation
#' of selection probabilities.
#' @param censoring censoring can get one of the three values "after" (default), "before", "none".
#' @param bs a logical value indicating whether to perform bootstap for estimation of
#' a confidence interval, and a SE. The default if FALSE.
#' If bs=FALSE there will not be provided any estimate of uncertainty.
#' @param nbs.rep number of bootstrap replications. The default is 200.
#' @param conf.int The confidence level for confidence intervals and hypotheses tests.
#' The default level is 0.95.
#' @param seed for the bootstrap
#'
#' @return  A list with components:
#' \tabular{llr}{
#' \code{time} \tab a vector of time points at which the survival curve has a step \tab   \cr
#' \code{surv} \tab the estimate of survival at time \code{time}+0. \tab  \cr
#' \code{SE} \tab the estimate of SE \tab   \cr
#' \code{CI.L} \tab the pointwise  of a lower limit of the confidence interval \tab  \cr
#' \code{CI.U} \tab the pointwise  of an upper limit of the confidence interval\tab  \cr}
#'
#' @examples
#' # under construction
#'
#' @references Vakulenko-Lagun, B. Qian, J. Chiou, S.-H. Wang, N. Betensky, R.A. 2020. Estimation under covariate-induced dependent truncation using inverse probability weighting. Submitted.
#'
#' @export
survfit.CW <- function(X, T, delta, Z, censoring="after" , bs=FALSE, nbs.rep=200)
{
  # Z is a matrix (n X p) or a data frame
  q <- dim(Z)
  n <- length(T)
  Z <- as.matrix(Z, q[1], q[2])
  # removing missing observations:
  na.i <- is.na(X) | is.na(T) | is.na(delta) | apply(is.na(Z), 1, any)
  if (sum(na.i)>0)
  {
    cat(sum(na.i), " missing observations were omitted.\n")
    X <- X[!na.i]
    T <- T[!na.i]
    delta <- delta[!na.i]
    Z <- Z[!na.i,]
  }

  if (censoring=="after")
    est.CW <- .get.CW.res.cens(X, T, delta, Z)
  if (censoring %in% c("before", "none"))
    est.CW <- .get.CW(X, T, delta, Z)

  # bootstrap:
  # estimate SE at these time points:
  x=est.CW$time
  se.bs=cil=ciu=NULL

  if (bs) {
    b.b <-  NULL
    for (b in 1:nbs.rep)
    {
      samp.b <- sample(n, size=n, replace=TRUE)
      if (censoring=="after")
        bs.est <- .get.CW.res.cens(X[samp.b], T[samp.b], delta[samp.b], Z[samp.b,])
      if (censoring %in% c("before", "none"))
        bs.est <- .get.CW(X[samp.b], T[samp.b], delta[samp.b], Z[samp.b,])
      res.bs <- sapply(x, .get.S, time=bs.est$time, surv=bs.est$surv, simplify=TRUE, USE.NAMES=FALSE)
      b.b <- rbind(b.b, res.bs)
    }
    se.bs <- apply(b.b, 2, sd, na.rm=TRUE)
    cil <- apply(b.b, 2, quantile, prob=0.025, na.rm=TRUE)
    ciu <- apply(b.b, 2, quantile, prob=0.975, na.rm=TRUE)
  }
  list(time=est.CW$time, surv=est.CW$surv, SE=se.bs, CI.L=cil, CI.U=ciu)
}


# inner function
.get.TVW.res.cens <- function(X, T, delta, Z)
{
  f.C <- survfit(Surv(X-T, 1-delta) ~ 1)
  S_c <- sapply(X[delta==1]-T[delta==1], .get.S, time=summary(f.C)$time, surv=summary(f.C)$surv, simplify=TRUE, USE.NAMES=FALSE)

  # using only uncensored observations estimate truncation distribution given Z
  nn <- sum(delta)
  del= rep(1, nn)
  m <- max(c(X[delta==1], T[delta==1]))
  srv.trunc <- Surv(time=m-X[delta==1], time2=m-T[delta==1], event=del)
  #   sol <- try(coxph(srv.trunc ~ Z[delta==1,], weights=1/S_c) ) # if Z is a matrix....
  sol <- try(coxph(srv.trunc ~ Z[delta==1], weights=1/S_c) )
  if (class(sol)!="coxph")
  {
    cat("coxph: Error occurred in 'get.CW.res.cens'.","\n")
  }
  alpha <- matrix(sol$coef, length(sol$coef), 1)

  bh <- basehaz(sol, centered=FALSE)
  H <- bh$hazard # ordered cumulative hazard
  K.S.NA <- exp(-H) # in the order of m-T[n], m-T[n-1], ...
  n1 <- length(bh$hazard)

  K.CDF.forw <- K.S.NA[n1:1] # this is only baseline pr-ties for z=0!!!!
  T.times <- m-bh$time[n1:1]

  # begin Jing's -------------------
  X.s1 <- X[delta==1]
  Z.s1 <- Z[delta==1]
  or <- order(X.s1)
  n <- length(X)
  n1 <- sum(delta) # the number of events (including duplicates)

  Xj <- matrix(X, n, n1, byrow = FALSE)
  Xi <- matrix(X.s1[or], n, n1, byrow = TRUE)
  Tj <- matrix(T, n, n1, byrow = FALSE)
  Zj <- matrix(Z, n, n1, byrow = FALSE)

  base.jing <- 1/(1-sapply(X.s1[or], .get.S, time=T.times, surv=1-K.CDF.forw, simplify=TRUE, USE.NAMES=FALSE))
  w.jing <- base.jing^exp(c(alpha)*Z.s1[or])

  K <-  matrix(base.jing, n, n1, byrow = TRUE)^exp(c(alpha)*Zj)
  at.risk <- colSums((Tj<=Xi)*(Xi<=Xj)*K ) # vector 1Xncol
  # S.jing <- cumprod(1-w.jing/at.risk)
  # end Jing's --------------------

  tied <- data.frame(time=X.s1[or], base=base.jing, numer=w.jing, at.risk=at.risk)
  agg <- aggregate(numer ~ time, data = tied, sum)
  at.risk.u <- aggregate(at.risk ~ time, data = tied, max)

  S.jing <- cumprod(1-agg$numer/at.risk.u$at.risk)

  list(time=agg$time, surv=S.jing)
}


# inner function for TVW, allows ties(!!): use this function:
# !!!!! only for one covariate - I have to extend it to multiple covariates
.get.my.TVW <- function(X, T, delta, Z)
{
  d.t <- rep(1, length(X))
  srv.trunc <- Surv(time=m-X, time2=m-T, event=d.t)
  sol <- try(coxph(srv.trunc ~ Z) )
  if (class(sol)!="coxph")
  {
    cat("coxph: Error occurred. rep=", j,"\n")
  }
  alpha <- sol$coef[1]

  bh <- basehaz(sol, centered=FALSE)
  H <- bh$hazard # ordered cumulative hazard
  K.S.NA <- exp(-H) # in the order of m-T[n], m-T[n-1], ...
  n1 <- length(bh$hazard)

  K.CDF.forw <- K.S.NA[n1:1] # this is only baseline pr-ties for z=0!!!!
  T.times <- m-bh$time[n1:1]

  # from residual treating ties:
  X.s1 <- X[delta==1]
  Z.s1 <- Z[delta==1]
  or <- order(X.s1)
  n <- length(X)
  n1 <- sum(delta) # the number of events (including duplicates)

  Xj <- matrix(X, n, n1, byrow = FALSE)
  Xi <- matrix(X.s1[or], n, n1, byrow = TRUE)
  Tj <- matrix(T, n, n1, byrow = FALSE)
  Zj <- matrix(Z, n, n1, byrow = FALSE)

  base.jing <- 1/(1-sapply(X.s1[or], .get.S, time=T.times, surv=1-K.CDF.forw, simplify=TRUE, USE.NAMES=FALSE))
  w.jing <- base.jing^exp(c(alpha)*Z.s1[or])

  K <-  matrix(base.jing, n, n1, byrow = TRUE)^exp(c(alpha)*Zj)
  at.risk <- colSums((Tj<=Xi)*(Xi<=Xj)*K ) # vector 1Xncol

  tied <- data.frame(time=X.s1[or], base=base.jing, numer=w.jing, at.risk=at.risk)
  agg <- aggregate(numer ~ time, data = tied, sum)
  at.risk.u <- aggregate(at.risk ~ time, data = tied, max)

  S.jing <- cumprod(1-agg$numer/at.risk.u$at.risk)

  list(time=agg$time, surv=S.jing)
}



#' Time-varying-weight estimation of survival probability
#'
#'
#' The estimator assumes the Cox proportional hazards model for truncation time given covariates, which is
#' used to obtain "time-varying weights" and reweigh the contribution of the observations at-risk accordingly.
#' For further details we refer to the paper of Vakulenko-Lagun et al. (2020).
#'
#'
#' @param X a vector of observed (possibly right-censored) lifetimes.
#' @param T a vector of truncation times.
#' @param delta a censoring indicator. It is assumed that 0 corresponds to censoring.
#' @param Z a matrix of covariates (common causes of \code{T} and \code{X}) that will be used for estimation
#' of selection probabilities. NOTE: so far implemented only for one covariate.
#' @param censoring censoring can get one of the three values "after" (default), "before", "none".
#' @param bs a logical value indicating whether to perform bootstap for estimation of
#' a confidence interval, and a SE. The default if FALSE.
#' If bs=FALSE there will not be provided any estimate of uncertainty.
#' @param nbs.rep number of bootstrap replications. The default is 200.
#' @param conf.int The confidence level for confidence intervals and hypotheses tests.
#' The default level is 0.95.
#' @param seed for the bootstrap
#'
#' @return  A list with components:
#' \tabular{llr}{
#' \code{time} \tab a vector of time points at which the survival curve has a step \tab   \cr
#' \code{surv} \tab the estimate of survival at time \code{time}+0. \tab  \cr
#' \code{SE} \tab the estimate of SE \tab   \cr
#' \code{CI.L} \tab the pointwise  of a lower limit of the confidence interval \tab  \cr
#' \code{CI.U} \tab the pointwise  of an upper limit of the confidence interval\tab  \cr}
#'
#' @examples
#' # under construction
#'
#' @references Vakulenko-Lagun, B. Qian, J. Chiou, S.-H. Wang, N. Betensky, R.A. 2020. Estimation under covariate-induced dependent truncation using inverse probability weighting. Submitted.
#'
#' @export
survfit.TVW <- function(X, T, delta, Z, censoring="after" , bs=FALSE, nbs.rep=200)
{
  # Z is a matrix (n X p) or a data frame
  q <- dim(Z)
  n <- length(T)
  Z <- as.matrix(Z, q[1], q[2])
  # removing missing observations:
  na.i <- is.na(X) | is.na(T) | is.na(delta) | apply(is.na(Z), 1, any)
  if (sum(na.i)>0)
  {
    cat(sum(na.i), " missing observations were omitted.\n")
    X <- X[!na.i]
    T <- T[!na.i]
    delta <- delta[!na.i]
    Z <- Z[!na.i,]
  }

  if (censoring=="after")
    est.TVW <- .get.TVW.res.cens(X, T, delta, Z)
  if (censoring %in% c("before", "none"))
    est.TVW <- .get.my.TVW(X, T, delta, Z)

  # bootstrap:
  # estimate SE at these time points:
  x=est.TVW$time
  se.bs=cil=ciu=NULL

  if (bs) {
    b.b <-  NULL
    for (b in 1:nbs.rep)
    {
      samp.b <- sample(n, size=n, replace=TRUE)
      if (censoring=="after")
        bs.est <- .get.TVW.res.cens(X[samp.b], T[samp.b], delta[samp.b], Z[samp.b,])
      if (censoring %in% c("before", "none"))
        bs.est <- .get.my.TVW(X[samp.b], T[samp.b], delta[samp.b], Z[samp.b,])
      res.bs <- sapply(x, .get.S, time=bs.est$time, surv=bs.est$surv, simplify=TRUE, USE.NAMES=FALSE)
      b.b <- rbind(b.b, res.bs)
    }

    se.bs <- apply(b.b, 2, sd, na.rm=TRUE)
    cil <- apply(b.b, 2, quantile, prob=0.025, na.rm=TRUE)
    ciu <- apply(b.b, 2, quantile, prob=0.975, na.rm=TRUE)
  }
  list(time=est.TVW$time, surv=est.TVW$surv, SE=se.bs, CI.L=cil, CI.U=ciu)
}


