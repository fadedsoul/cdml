
################################################################################
#             model: continuous treatment effect / IV                          #
#  response.y.: polynom case, with normal noise, covariates, treatment         #
#                    g.method <- c("lasso")                                    #
#          gps.method <- c( "boosting&normal", "linear&boxcox")                #
################################################################################

################ FUNCTION ###########################

#' Simulation of using lasso for linear and nonlinear structure model
#'
#' @param model choose from "IV" and "CTE"
#'
#' @return files with certain name tags, e.g. "linear.lasso.cte.Rdata"
#' @export
application.lasso <- function(model) {

  if (model == "IV") {
    file1 <- "linear.lasso.iv.Rdata"
    file2 <- "polynom.lasso.iv.Rdata"
    file3 <- "polynom2.lasso.iv.Rdata"
    file4 <- "polynom3.lasso.iv.Rdata"
  }
  if(model == "CTE") {
    file1 <- "linear.lasso.cte.Rdata"
    file2 <- "polynom.lasso.cte.Rdata"
    file3 <- "polynom2.lasso.cte.Rdata"
    file4 <- "polynom3.lasso.cte.Rdata"
  }

  running.simulation(
    model <- model,
    simu <- 40,
    samples <- c(200, 500, 1000, 2000),
    g.method <- c("lasso"),
    gps.method <-
      c(#"series",
        "rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "linear",
    sd = 8,
    file <- file1
  )

  running.simulation(
    model <- model,
    simu <- 40,
    samples <- c(200, 500, 1000, 2000),
    g.method <- c("lasso"),
    gps.method <-
      c(#"series",
        "rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd = 8,
    file <- file2 # the file saves data
  )

  running.simulation(
    model <- model,
    simu <- 40,
    samples <- c(200, 500, 1000, 2000),
    g.method <- c("lasso"),
    gps.method <-
      c(#"series",
        "rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom2",
    sd = 8,
    file <- file3 # the file saves data
  )

  running.simulation(
    model <- model,
    simu <- 40,
    samples <- c(200, 500, 1000, 2000),
    g.method <- c("lasso"),
    gps.method <-
      c(#"series",
        "rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom3",
    sd = 8,
    file <- file4 # the file saves data
  )
}

