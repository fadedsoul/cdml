
################################################################################
#                                                                              #
#                  Purpose:         true data generation                       #
#                                                                              #
################################################################################

#' True Data Generation
#'
#' @param numOfCovariates dimension of covariates
#' @param distributionForCovariates default = "normal", choose from normal, student, uniform distributions
#' @param distributionForNoise default = "normal", choose from normal, student, uniform distributions
#' @param distributionForTreatment default = "normal", choose from normal, student, uniform distributions
#' @param parametricCurveOption default = "linear", choose from "linear", "polynom", "polynom2", "polynom3", "mixture"
#' @param trimLowerBound the lower bound of trimming treatment vector
#' @param trimUpperBound the upper bound of trimming treatment vector
#' @param nrep default = 1000000, numerically calculate the true density by setting equi-distant nrep sample points for approximation
#'
#' @return a list of true response value "true.y" evaluated at "grid", with treatment density "density".
#' @export
#'
true.data.generation <- function(numOfCovariates = 2,
                                 distributionForCovariates = "normal",
                                 distributionForNoise = "normal",
                                 distributionForTreatment = "normal",
                                 parametricCurveOption = "polynom",
                                 trimLowerBound = -4,
                                 trimUpperBound = 4,
                                 nrep = 1000000)
{
  xx <-
    covariates.generation(
      numOfCovariates = numOfCovariates,
      numOfSamples = nrep,
      distributionOption = distributionForCovariates
    )
  tt <-
    treatment.generation(
      covariates = xx,
      numOfSamples = nrep,
      distributionOption = distributionForTreatment
    )

  data2 <- data.frame(xx, tt)
  data2 <-
    data2[(data2$tt > trimLowerBound & trimUpperBound > data2$tt),]
  row.names(data2) <- 1:nrow(data2)

  n = ceiling(100 * abs(trimLowerBound - trimUpperBound))

  true.density <- density(data2$tt,
                          from = trimLowerBound,
                          to = trimUpperBound,
                          n = n)

  if (parametricCurveOption == "linear") {
    tt <- seq(from = trimLowerBound,
              to = trimUpperBound,
              length.out = n)
    true.y <- tt + mean(rowSums(xx))
  }
  if (parametricCurveOption == "polynom") {
    tt <- seq(from = trimLowerBound,
              to = trimUpperBound,
              length.out = n)
    true.y <- (1 / 2) * tt ^ 2 + tt
  }
  if (parametricCurveOption == "polynom2") {
    tt <- seq(from = trimLowerBound,
              to = trimUpperBound,
              length.out = n)
    true.y <-
      (1 / 4) * tt ^ 3+(1 / 2) * tt ^ 2 + tt
  }
  if (parametricCurveOption == "polynom3") {
    tt <- seq(from = trimLowerBound,
              to = trimUpperBound,
              length.out = n)
    true.y <- (1 / 2) * tt ^ 2 + tt
  }

  newList <- data.frame("grid" = true.density$x,
                        "density" = true.density$y,
                        "true.y" = true.y)

  return(newList)
}
