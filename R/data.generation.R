


################################################################################
#                                                                              #
#                Wrap Up the function data.generation                          #
#                     to run in parallel on many CPUs                          #
#                                                                              #
################################################################################

#' Data Generation wrap-up (short-cut name: dg.w)
#'
#' @param simu number of times one would run, default = 100.
#' @param samples number of samples, default = 1000
#' @param covariates dimension of covariates
#' @param model choose from "CTE", "IV"
#' @param trimLowerBound the lower bound of treatment trimming t, default = -4.
#' @param trimUpperBound the lower bound of treatment trimming t, default = 4.
#' @param parametricCurveOption choose from "linear", "polynom", "polynom2", "polynom3", "mixture".
#' @param ncores number of cores we use
#' @param sd  \eqn{t = \sum x_i + \epsilon}, the standard error of \eqn{\epsilon}, choose from 1, 2, 3, 5, 8, 10, 15
#' @param noise.sd the standard error of response generation assuming the response model is based on a normal distribution
#'
#' @return list that contains two parts:
#'                 $data contains generated data,
#'                 $true contains the true underlying response and the density at grid on treatment
#' @export


data.generation.wrap <-
  function(simu = 100,
           samples = 1000,
           covariates = 2,
           model,
           trimLowerBound = -4,
           trimUpperBound = 4,
           parametricCurveOption = "polynom",
           ncores = ncores,
           sd = 8,
           noise.sd = 1)
  {
     data.list <- parallel::mcmapply(
      numOfSamples <- rep(samples, times = simu),
      FUN = data.generation,

      MoreArgs = list(
        noise.sd = noise.sd,
        sd=sd,
        numOfCovariates = covariates,
        model = model,
        trimLowerBound = trimLowerBound,
        trimUpperBound = trimUpperBound,
        parametricCurveOption = parametricCurveOption,
        numOfPreSamples = NULL
      ),
      SIMPLIFY = FALSE,
      mc.cores = ncores
    )

    return(data.list)
  }


dg.w <- data.generation.wrap

################################################################################
#                                                                              #
#                  Purpose:           data generation                          #
#                    (Merging dataframe and Trimming)                          #
#                                                                              #
################################################################################

#'  Data Generation
#'
#' @param model choose from "CTE", "IV"
#' @param numOfSamples numer of samples, default = 1000.
#' @param numOfPreSamples number of samples before sample trimming. As default, we set this as NULL
#' @param numOfCovariates dimension of covariates, default = 2.
#' @param distributionForCovariates choose from normal, uniform, student distributions.
#' @param distributionForIVs choose from normal, uniform, student distributions.
#' @param distributionForNoise choose from normal, uniform, student distributions.
#' @param distributionForTreatment choose from normal, uniform, student distributions.
#' @param parametricCurveOption choose from "linear", "polynom", "polynom2", "polynom3", "mixture".
#' @param trimLowerBound the lower bound of treatment trimming t.
#' @param trimUpperBound the upper bound of treatment trimming t.
#' @param sd  \eqn{t = \sum x_i + \epsilon}, the standard error of \eqn{\epsilon}
#' @param noise.sd the standard error of response generation assuming the response model is based on a normal distribution
#'
#' @return frame of generated data.
#' @export
#'
#' @examples data.generation(model = "CTE", numOfSamples = 1000,
#' numOfCovariates = 2,
#' distributionForCovariates = "normal",
#' distributionForNoise = "normal",
#' distributionForTreatment = "normal",
#' parametricCurveOption = "linear",
#' trimLowerBound = -10,
#' trimUpperBound = 10,
#' sd = 8)
data.generation <- function(model = "CTE",
                            numOfSamples = 1000,
                            numOfPreSamples = NULL,
                            numOfCovariates = 2,
                            distributionForCovariates = "normal",
                            distributionForIVs = "normal",
                            distributionForNoise = "normal",
                            distributionForTreatment = "normal",
                            parametricCurveOption = "linear",
                            trimLowerBound = -10,
                            trimUpperBound = 10,
                            sd,
                            noise.sd = 1)
{
  if ((!is.null(numOfPreSamples)) && (!is.null(numOfSamples))) {
    stop("Please enter either numOfSamples or numOfPreSamples in function data.generation!")
  }


  if (!is.null(numOfSamples)) {
    numOfPreSamples <- ceiling(10 * numOfSamples)
    #numOfPreSamples <- ceiling(100 * numOfSamples)
  }

  ### if the model we choose here is continuous treatment effect model, as default
  if(model == "CTE"){
    x <-
      covariates.generation(
        numOfCovariates = numOfCovariates,
        numOfSamples = numOfPreSamples,
        distributionOption = distributionForCovariates
      )
    u <-
      noise.generation(numOfSamples = numOfPreSamples, distributionOption = distributionForNoise, noise.sd = noise.sd)
    t <-
      treatment.generation(
        covariates = x,
        numOfSamples = numOfPreSamples,
        distributionOption = distributionForTreatment,
        sd=sd
      )
    y <- response.generation(
      covariates = x,
      treatment = t,
      noise = u,
      numOfSamples = numOfPreSamples,
      parametricCurveOption = parametricCurveOption
    )

    data <- data.frame(t, x, y)
    data <-
      data[(data$t > trimLowerBound & trimUpperBound > data$t),]
    row.names(data) <- 1:nrow(data)
  }

  ### if the model we choose here is intrumental variable model
  if(model == "IV"){
    x <-
      covariates.generation(
        numOfCovariates = numOfCovariates,
        numOfSamples = numOfPreSamples,
        distributionOption = distributionForCovariates
      )

    z <- instrument.generation(
      covariates = x,
      numOfSamples = numOfPreSamples,
      distributionOption = distributionForIVs
    )

    u <-
      noise.generation(numOfSamples = numOfPreSamples, distributionOption = distributionForNoise)
    t <-
      treatment.generation(
        covariates = x,
        numOfSamples = numOfPreSamples,
        distributionOption = distributionForTreatment      )
    y <- response.generation(
      covariates = x,
      treatment = t,
      noise = u,
      numOfSamples = numOfPreSamples,
      parametricCurveOption = parametricCurveOption
    )



    data <- data.frame(t, z, x, y)
    data <-
      data[(data$t > trimLowerBound & trimUpperBound > data$t),]
    row.names(data) <- 1:nrow(data)

  }
  if (!is.null(numOfSamples)) {
    if (dim(data)[1] >= numOfSamples) {
      data <- data[1:numOfSamples, ]
    } else{
      stop("The number of the samples is less than required!")
    }
  }

  return(data)
}
