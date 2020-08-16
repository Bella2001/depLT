#' Package for Estimation of Survival Function under Dependent Left Truncation
#'
#' The \pkg{depLT} package accompanies the paper of Vakulenko-Lagun et al. (2020).
#' It is designed for analysis of left-truncated data where truncation time and lifetime are dependent.
#' We assume that the dependence is induced by covariates, or common causes of the lifetime and truncation time,
#' and that all these common causes are measured.
#' The package implements two methods, a case-weights (CW) estimator and a time-varying-weights (TVW) estimator, and
#' two types of censoring, residual censoring that can happen only after sampling and an original-time-scale censoring,
#' that can occur before sampling into the study.
#'
#' The \pkg{depLT} package provides two main functions,
#' \code{\link{survfit.CW}}
#'  and \code{\link{survfit.TVW}}, and
#'  functions \code{\link{positivity.binary}} and \code{\link{positivity.cont}} for diagnostics of practical, estimator-specific violation of positivity.

#'
#'
#' @references Vakulenko-Lagun, B. Qian, J. Chiou, S.-H. Wang, N. Betensky, R.A. 2020. Estimation under covariate-induced dependent truncation using inverse probability weighting. Submitted.
#'
#' @import survival
#' @import tidyverse
#' @import inline

#'
#'
#' @importFrom stats pnorm qnorm quantile sd var
#' @useDynLib depLT, .registration = TRUE
#'
#' @docType package
#' @name depLT
NULL
