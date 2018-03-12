################################################################################
#                                                                              #
#                Wrap Up the function simulation                               #
#                     to run in parallel on many CPUs                          #
#                                                                              #
################################################################################

#' Simulation wrap-up (short-cut name: sim.w)
#'
#' @param kfold number of folds for sample splitting
#' @param data data
#' @param model model varies from "IV", "CTE"
#' @param ncores number of cores we use
#' @param method estimation method we choose from double machine learning method ("CDML")
#'            and simple regression method ("SR"), Hirano & Imbens method ("HI")
#' @param g.method method for regression estimation
#' @param gps.method method for generalized propensity score estimation
#' @param trimLowerBound.t trimming lower bound of treatment vector
#' @param trimUpperBound.t trimming upper bound of treatment vector
#'
#' @return a list of returned objects in function simulation.
#' @export
#'
simulation.wrap <-
  function(kfold,
           data,
           model,
           ncores = ncores,
           method,
           g.method = "rf",
           gps.method = "series",
           trimLowerBound.t = -4,
           trimUpperBound.t = 4,
           detoured = FALSE)
  {
    result.list <- mcmapply(
      data <- data,
      FUN = simulation,
      MoreArgs = list(
        kfold = kfold,
        trimLowerBound.pi = 0.0001,
        trimLowerBound.pibar = 0.0001,
        trimUpperBound.pi = 1,
        trimUpperBound.pibar = 1,
        model = model,
        method = method,
        g.method = g.method,
        gps.method = gps.method,
        trimLowerBound.t = trimLowerBound.t,
        trimUpperBound.t = trimUpperBound.t,
        variation = "whole",
        verbose = FALSE,
        detoured = detoured
      ),
      SIMPLIFY = FALSE,
      mc.cores = ncores
    )

    ## for keeping track in the big-simulation
    print("fire!!")

    return(result.list)
  }
sim.w <- simulation.wrap
