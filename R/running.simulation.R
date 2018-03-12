
################################################################################
#                                                                              #
#                         Simulation Running                                   #
#                                                                              #
################################################################################

#' Running Simulation
#'
#' @param model model varies from "IV", "CTE"
#' @param simu number of simulations one would run, default = 100.
#' @param samples number of samples
#' @param g.method a vector of method for regression estimation
#' @param gps.method a vector of method for generalized propensity score estimation
#' @param trimming trimming of treatment vector, default = c(-4, 4)
#' @param cov a vector of dimensions of covariates
#' @param method a vector of estimation methods we choose from double machine learning method ("CDML")
#'            and simple regression method ("SR"), Hirano & Imbens method ("HI")
#' @param fold number of folds for sample splitting
#' @param responseCurve choose from "linear", "polynom", "polynom2", "polynom3", "mixture".
#' @param file file that we save our return into
#' @param sd \eqn{t = \sum x_i + \epsilon}, the standard error of \eqn{\epsilon}, choose from 1, 2, 3, 5, 8, 10, 15
#'
#' @return no returns
#' @export
#'
#'
#'
#'
#' @import parallel
#' @import randomForest
#' @import KernSmooth
#' @import stats
#' @import nnet
#' @import rpart
#' @import earth
#' @import glmnet
#' @import EnvStats
#' @import gbm
#' @import polycor
#' @import np
#' @import utils
#' @import quantregForest
#'
#' @examples running.simulation(
#' ########### enter all parameters needed for running simulations ##########
#' model <- "CTE",
#' simu <-8,
#' samples <- c(200),
#' g.method <- c("lasso"),
#' gps.method <-
#'   c(#"series",
#'    "rf&normal" ,
#'    "linear&boxcox"),
#' trimming <- c(-4, 4),
#'cov <- c(5),
#' method <- c("SR", "CDML"),
#' fold <- c(1, 2, 3),
#' responseCurve <- "polynom3",
#' file <- "demo.Rdata", # the file saves data
#' sd =8)
running.simulation <-
  function(model = "CTE",
           simu = 40,
           samples = c(100, 200, 500),
           g.method = c("nnet", "rf", "MARS"),
           # i deleted lasso for simplification, may discuss that case later.
           gps.method =
             c(#"nnet&normal",
               #"rf&normal",
               #"nnet&boxcox",
               #"rf&boxcox",
               #"linear&normal",
               "series",
               "boosting&normal",
               "linear&boxcox"),
           trimming = c(-4, 4),
           cov = c(2, 5, 10),
           method = c("SR", "HI", "CDML"),
           fold = c(1, 2, 3, 5),
           responseCurve = "polynom",
           file =
             "r_scratch/cte100and200and500.Rdata",
           sd = 8)
  {

    set.seed(101) # for reproducibility

    ncores <-
      min(detectCores(all.tests = FALSE, logical = TRUE), 44) # CPUs that we use

    simu.setups.inbetween1 <- NULL
    simu.setups.inbetween2 <- NULL
    simu.setups.inbetween3 <- NULL
    ############################ Set-Ups #########################################
    if ("SR" %in% method) {
     simu.setups.inbetween1 <-
      expand.grid(
        list(
          "g.method" =  g.method,
          "gps.method" = "linear&normal",
          "method" = "SR",
          "fold" = 1
        )
      )
    }
    if ("CDML" %in% method) {
    simu.setups.inbetween2 <-
      expand.grid(
        list(
          "g.method" =  g.method,
          "gps.method" = gps.method,
          "method" = "CDML",
          "fold" = fold
        )
      )
    }

    if ("HI" %in% method) {
      simu.setups.inbetween3 <-
        expand.grid(list(
          "g.method" =  "nnet",
          "gps.method" = gps.method,
          "method" = "HI",
          "fold" = 1
        ))
    }
    simu.setups.inbetween <-
      rbind(simu.setups.inbetween1,
            simu.setups.inbetween2,
            simu.setups.inbetween3)
    simu.setups.inbetween <-
      simu.setups.inbetween[which(((simu.setups.inbetween$fold != 1) &
                                     (simu.setups.inbetween$method != "CDML")
      ) == FALSE),]
    row.names(simu.setups.inbetween) <-
      1:nrow(simu.setups.inbetween)

    data.setups.inbetween <-
      expand.grid(list("samples" = samples,
                       "covariates" = cov))


    simu.setups <- simu.setups.inbetween
    data.setups <- data.setups.inbetween

    #########################################################################

    print("fire, taget number:")
    cat("\n")
    print((dim(data.setups)[1]) * (dim(simu.setups)[1]))

    #########################################################################

    data.list.list <-
      mapply(
        simu = simu,
        samples = data.setups$samples,
        covariates = data.setups$covariates,
        FUN = data.generation.wrap,
        MoreArgs = list(
          trimLowerBound = trimming[1],
          trimUpperBound = trimming[2],
          parametricCurveOption = responseCurve,
          model = model,
          ncores = ncores,
          sd = sd
        ),
        SIMPLIFY = FALSE
      )
    #################################################
    true.data <-
      true.data.generation(
        parametricCurveOption = responseCurve,
        trimLowerBound = trimming[1],
        trimUpperBound = trimming[2]
      )
    ##################################################
    res.list.list <- mapply(
      data = rep(data.list.list, each = length(simu.setups$method)),
      method = rep(simu.setups$method, times = length(data.list.list)),
      kfold = rep(simu.setups$fold, times = length(data.list.list)),
      g.method = rep(simu.setups$g.method, times = length(data.list.list)),
      gps.method = rep(simu.setups$gps.method, times = length(data.list.list)),
      FUN = simulation.wrap,
      MoreArgs = list(
        trimLowerBound.t = trimming[1],
        trimUpperBound.t = trimming[2],
        model = model,
        ncores = ncores
      ),
      SIMPLIFY = FALSE
    )

    ##############################################
    save(true.data,
         simu.setups,
         data.setups,
         res.list.list,
         file = file)
    #################################################
    error.list.list <- mapply(
      result.list = res.list.list,
      FUN = err.w,
      MoreArgs = list(true.data = true.data, ncores = ncores),
      SIMPLIFY = FALSE
    )
    # ################################################

    save(true.data,
         simu.setups,
         data.setups,
         res.list.list,
         error.list.list,
         file = file)

  }

if(FALSE) {
    ## Alternative  to  save()  and  load()
    saveRDS(list(true.data =true.data,
                 simu.setups = simu.setups,
                 data.setups = data.setups,
                 res.list.list = res.list.list),
            file = file) ## << convention use file "<something>.rds"

    ## use as
    sim.123 <- readRDS(file)
}
