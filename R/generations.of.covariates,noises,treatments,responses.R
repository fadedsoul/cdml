################################################################################
#                                                                              #
#   Purpose: Generation of Covariates, Noises, Treatments and Responses        #
#                                                                              #
################################################################################

#' Generation of the covariates dataset
#'
#' @param distributionOption choose from normal, uniform, student distributions
#' @param numOfCovariates dimension of covariates
#' @param numOfSamples number of samples
#'
#' @return A matrix the columns of which represent different covariates, the rows of which represents different observations.
#'
covariates.generation <-
  function(distributionOption = "normal",
           numOfCovariates = 2,
           numOfSamples)
  {
    check.distribution.name <- function(name) {
      return(name == "normal" |
               name == "uniform" |
               name == "student" | name == "mixed")
    }

    if (!check.distribution.name(distributionOption)) {
      stop("This distribution that you have input for the covariates is not defined!")
    }
    ## here we generate the covariates in a way that covariates are independent with each other
    if (distributionOption == "normal") {
      xmat <-
        matrix(
          rnorm(
            numOfSamples * numOfCovariates,
            mean = 0,
            sd = 1
          ),
          nrow = numOfSamples,
          ncol = numOfCovariates
        )
    }

    if (distributionOption == "uniform") {
      xmat <-
        matrix(
          runif(
            numOfSamples * numOfCovariates,
            min = -1,
            max = 1
          ),
          nrow = numOfSamples,
          ncol = numOfCovariates
        )
    }

    if (distributionOption == "student") {
      xmat <-
        matrix(rt(numOfSamples * numOfCovariates,
                  df = 3),
               nrow = numOfSamples,
               ncol = numOfCovariates)
    }

    if (distributionOption == "mixed") {
      xmat <- matrix(nrow = numOfSamples,
                     ncol = numOfCovariates)
      for (i in 1:numOfCovariates) {
        if (i %% 3 == 0) {
          xmat[, i] = rnorm(numOfSamples, mean = 0, sd = 1)
        } else if (i %% 3 == 1) {
          xmat[, i] = runif(numOfSamples, min = -1, max = 1)
        } else if (i %% 3 == 2) {
          xmat[, i] = rt(numOfSamples, df = 3)
        }
      }
    }

    return(xmat)
  }

#' Generation of treatment dataset
#'
#' @param distributionOption choose from normal, uniform, student distributions
#' @param numOfSamples number of samples
#' @param sd \eqn{t = \sum x_i + \epsilon}, the standard error of \eqn{\epsilon}
#' @param covariates the object that returned from function covariates.generation
#'
#' @return A vector stands for treatment
#'
treatment.generation <-
  function(distributionOption = "normal",
           numOfSamples = 1000,
           covariates,
           sd=8)
  {
    check.distribution.name <- function(name) {
      return(
        name == "normal" |
          name == "uniform" |
          name == "student" |
          name == "mixed"
      )
    }

    if (!check.distribution.name(distributionOption)) {
      stop("This distribution that you have input for the treatment is not defined!")
    }

    if (distributionOption == "normal") {
      tmat <-
        apply(covariates, 1, function(x)
          rnorm(1, mean = sum(x), sd = sd))
    }

    if (distributionOption == "uniform") {
      tmat <-
        apply(covariates, 1, function(x)
          runif(1, min = -0.5, max = 0.5))
    }

    if (distributionOption == "student") {
      tmat <-
        apply(covariates, 1, function(x)
          rt(1, df = 3))
    }

    return(tmat)
  }

#' Generation of noise dataset
#'
#' @param distributionOption choose from normal, uniform, student distributions
#' @param numOfSamples number of samples
#' @param noise.sd if our noise is normal distributed, then we use this parameter to tune the noise. Default value equals 1
#'
#' @return A vector stands for noise
#'
noise.generation <-
  function(distributionOption = "normal",
           numOfSamples = 1000,
           noise.sd = 1)
  {
    check.distribution.name <- function(name) {
      return(name == "normal" |
               name == "uniform" |
               name == "student")
    }
    if (!check.distribution.name(distributionOption)) {
      stop("This distribution that you have input for noise term is not defined!")
    }

    nmat <- vector('numeric', numOfSamples)

    if (distributionOption == "normal") {
      nmat <- rnorm(numOfSamples, mean = 0, sd = noise.sd)
    }

    if (distributionOption == "uniform") {
      nmat <- runif(numOfSamples, min = 0, max = 1)
    }

    if (distributionOption == "student") {
      nmat <- rt(numOfSamples, df = 3)
    }

    return(nmat)
  }


#' Generation of response dataset
#'
#' @param noise the noise vector generated by function noise.generation(), term u in the model
#' @param covariates the covariates matrix generated by function covariates.generation(), term x in the model
#' @param treatment the treatment vector generated by function treatment.generation(), term t in the model
#' @param parametricCurveOption g(x,t) with different smoothness chosen from "linear", "polynom", "polynom2", "polynom3", "mixture"
#' @param numOfSamples number of samples
#'
#' @return a vector of responses
#'
response.generation <-
  function(noise,
           covariates,
           treatment,
           parametricCurveOption = "linear",
           numOfSamples = 1000)
  {
    check.parametricCurve.name <- function(name) {
      return(
        name == "linear" |
          name == "polynom" |
          name == "polynom2" |
          name == "polynom3" | name == "mixture"
      )
    }
    if (!check.parametricCurve.name(parametricCurveOption)) {
      stop("This parametric model that you have input for the response term is not defined!")
    }

    rmat <- vector('numeric', numOfSamples)

    if (parametricCurveOption == "linear") {
      rmat <- treatment + rowSums(covariates) + noise
    }

    if (parametricCurveOption == "polynom") {
      rmat <-
        (1 / 2) * treatment ^ 2 + treatment + rowSums(covariates) + noise
    }

    if (parametricCurveOption == "polynom2") {
      rmat <-
        (1 / 4) * treatment ^ 3 + (1 / 2) * treatment ^ 2 + treatment + 3 *
        rowSums(covariates)  + noise
    }

    if (parametricCurveOption == "polynom3") {
      rmat <-
        (1 / 2) * treatment ^ 2 + treatment + rowSums(covariates) - 3 *
        treatment * covariates[, 1] + 2 * treatment * covariates[, 2] + noise
    }

    if (parametricCurveOption == "mixture") {
      ### TBA
    }

    return(rmat)
  }

#' Generation of instrument variable dataset
#'
#' @param distributionOption choose from normal, uniform, student distributions
#' @param covariates dimension of covariates
#' @param numOfSamples number of samples
#'
#' @return a vector of instrumental variable
#'
instrument.generation <-
  function(distributionOption = "normal",
           covariates,
           numOfSamples)
  {
    check.distribution.name <- function(name) {
      return(name == "normal" |
               name == "uniform" |
               name == "student" | name == "mixed")
    }

    if (!check.distribution.name(distributionOption)) {
      stop("This distribution that you have input for the treatment is not defined!")
    }

    if (distributionOption == "normal") {
      zmat <-
        apply(covariates, 1, function(x)
          rnorm(1, mean = sum(x), sd = 8))
    }

    if (distributionOption == "uniform") {
      zmat <-
        apply(covariates, 1, function(x)
          runif(1, min = -0.5, max = 0.5))
    }

    if (distributionOption == "student") {
      zmat <-
        apply(covariates, 1, function(x)
          rt(1, df = 3))
    }

    return(zmat)
  }
