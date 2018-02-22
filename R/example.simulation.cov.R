
################################################################################
#                         model: IV/CTE model                                  #
#  response.y.: polynom case, with normal noise, covariates, treatment         #
#          g.method <- c("nnet", "rf", "MARS")                                 #
#           gps.method <- c("boosting&normal", "linear&boxcox")                #
################################################################################


################ FUNCTION ###########################
#' Simulation of applying different dimensions of covariates
#'
#' @param model choose from "IV" and "CTE"
#'
#' @return files with certain name tags, e.g. "linear.cov.cte.Rdata"
#' @export
#'
application.covariates <- function(model) {


  if (model == "IV") {
    file1 <- "linear.cov.iv.Rdata"
    file2 <- "polynom.cov.iv.Rdata"
    file3 <- "polynom2.cov.iv.Rdata"
    file4 <- "polynom3.cov.iv.Rdata"
  }
  if(model == "CTE") {
    file1 <- "linear.cov.cte.Rdata"
    file2 <- "polynom.cov.cte.Rdata"
    file3 <- "polynom2.cov.cte.Rdata"
    file4 <- "polynom3.cov.cte.Rdata"
  }

 running.simulation(
   model <- model,
   simu <- 40,
   samples <- c(200, 500,1000,2000),
   g.method <- c("rf", "nnet", "trees"),
   # i deleted lasso for simplification, may discuss that case later.
   gps.method <-
     c(#"series",
       "rf&normal" ,
       #"boosting&normal"
       "linear&boxcox"),
   trimming <- c(-4, 4),
   cov <- c(2, 3, 4, 5, 7, 10, 15, 20),
   method <- c("SR", "CDML", "HI"),
   fold <- c(1, 2, 3, 5),
   responseCurve <- "linear",
   sd = 8,
   file <- file1 # the file saves data
 )

running.simulation(
  model <- model,
  simu <- 40,
  samples <- c(200, 500,1000,2000),
  g.method <- c("rf", "nnet", "trees"),
  # i deleted lasso for simplification, may discuss that case later.
  gps.method <-
    c(#"series",
      "rf&normal" ,
      #"boosting&normal"
      "linear&boxcox"),
  trimming <- c(-4, 4),
  cov <- c(2, 3, 4, 5, 7, 10, 15, 20),
  method <- c("SR", "CDML", "HI"),
  fold <- c(1, 2, 3, 5),
  responseCurve <- "polynom",
  sd = 8,
  file <- file2 # the file saves data
)

running.simulation(
  model <- model,
  simu <- 40,
  samples <- c(200, 500,1000,2000),
  g.method <- c("rf", "nnet", "trees"),
  # i deleted lasso for simplification, may discuss that case later.
  gps.method <-
    c(#"series",
      "rf&normal" ,
      #"boosting&normal"
      "linear&boxcox"),
  trimming <- c(-4, 4),
  cov <- c(2, 3, 4, 5, 7, 10, 15, 20),
  method <- c("SR", "CDML", "HI"),
  fold <- c(1, 2, 3, 5),
  responseCurve <- "polynom2",
  sd = 8,
  file <- file3 # the file saves data
)

running.simulation(
  model <- model,
  simu <- 40,
  samples <- c(200, 500,1000,2000),
  g.method <- c("rf", "nnet", "trees"),
  # i deleted lasso for simplification, may discuss that case later.
  gps.method <-
    c(#"series",
      "rf&normal" ,
      #"boosting&normal"
      "linear&boxcox"),
  trimming <- c(-4, 4),
  cov <- c(2, 3, 4, 5, 7, 10, 15, 20),
  method <- c("SR", "CDML", "HI"),
  fold <- c(1, 2, 3, 5),
  responseCurve <- "polynom3",
  sd = 8,
  file <- file4 # the file saves data
)
}

##########################################################
#                                                        #
#    Remark: the plotting function is saved locally !    #
#          Ask if needed.                                #
##########################################################
