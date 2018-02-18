################################################################################
#                                                                              #
#                     Purpose: error evaluation                                #
#                                                                              #
################################################################################

#' l^1 error evaluation
#'
#' @param result result from simulation function
#' @param true.data result from true.data.generation function
#'
#' @return  L1 error
error1.evaluation <- function(result, true.data)
{
  ## L1 error
  error.1 <-
    sum(true.data$density * abs((true.data$true.y - result$yhat)) ^ 1)

  return(error.1)
}

#' l^2 error evaluation
#'
#' @param result result from simulation function
#' @param true.data result from true.data.generation function
#'
#' @return  L2 error
error2.evaluation <- function(result, true.data)
{
  ## L2 error
  error.2 <-
    sum(true.data$density * abs((true.data$true.y - result$yhat)) ^ 2)

  return(error.2)
}


################################################################################
#                                                                              #
#            Bias evaluation procedure given res.list.list                     #
#                                                                              #
################################################################################

#' Bias Evaluation
#'
#' @param result.list  result list from simulation.wrap function
#' @param true.data result from true.data.generation function
#'
#' @return bias as a numeric
bias.evaluation <- function(result.list, true.data) {
  yhat <- 0
  for (i in 1:length(result.list)) {
    yhat <- 1 / length(result.list) * (result.list[[i]]$yhat) + yhat
  }
  bias <-
    sum(true.data$density * abs((true.data$true.y - yhat)))

  return(list("bias" = bias,
              "record" = result.list[[1]]$record))
}

#' Bias Evaluation Wrap-up
#'
#' @param result.list.list  list of result lists from running.simulation function
#' @param true.data result from true.data.generation function
#'
#' @return a list of bias lists
bias.evaluation.wrap <- function(res.list.list, true.data){
  bias.list.list <- mapply(
    res.list.list,
    FUN = bias.evaluation,
    MoreArgs = list(true.data = true.data),
    SIMPLIFY = FALSE
  )
  return(bias.list.list)
}
################################################################################
#                                                                              #
#                     Purpose: error evaluation wrap                           #
#                                                                              #
################################################################################
#' Error Evaluation Wrap-up
#'
#' @param result.list  result list from simulation.wrap function
#' @param true.data result from true.data.generation function
#'
#' @return a list of error lists where "error1" stands for l^1 error, "error2" stands for "l^2" error
error.evaluation.wrap <-
  err.w <- function(true.data, result.list, ncores)
    ## result.list: result.list from simulation wrap function
    ## true.data: result from true.data.generation function

    ## return value $error1, $error2 is a vector!
  {
    error1 <- mcmapply(
      result = result.list,
      FUN = error1.evaluation,
      MoreArgs = list(true.data = true.data),
      mc.cores = ncores
    )

    error2 <- mcmapply(
      result = result.list,
      FUN = error2.evaluation,
      MoreArgs = list(true.data = true.data),
      mc.cores = ncores
    )

    return(list(
      "error1" = error1,
      "error2" = error2,
      "record" = result.list[[1]]$record
    ))
  }
