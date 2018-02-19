
################################################################################
#                         model: CTE model                                     #
#  response.y.: polynom case, with normal noise, covariates, treatment         #
#               g.method <- c("nnet", "rf")                                    #
#           gps.method <- c("boosting&normal", "linear&boxcox")                #
################################################################################


################################ FUNCTION #################################
#' Simulation of using different standard errors in treatment generation
#'
#'
#' @return files with certain name tags, e.g. "polynom.sd1.Rdata"
#' @export
#'
application.sd <- function() {

## we fix responseCurve always here as "polynom"
running.simulation(
  model <- "CTE",
  simu <- 40,
  samples <- c(300, 500, 750, 1000, 1500),
  g.method <- c("rf", "nnet"),
  gps.method <-
    c("rf&normal" ,
      "linear&boxcox"),
  trimming <- c(-4, 4),
  cov <- c(5),
  method <- c("SR", "CDML", "HI"),
  fold <- c(1, 2, 3, 5),
  responseCurve <- "polynom",
  sd <- 1,
  file <- "polynom.sd1.Rdata" # the file saves data
)

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 2,
    file <- "polynom.sd2.Rdata" # the file saves data
  )

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 3,
    file <- "polynom.sd3.Rdata" # the file saves data
  )

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 5,
    file <- "polynom.sd5.Rdata" # the file saves data
  )

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 8,
    file <- "polynom.sd8.Rdata" # the file saves data
  )

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 10,
    file <- "polynom.sd10.Rdata" # the file saves data
  )

  running.simulation(
    model <- "CTE",
    simu <- 40,
    samples <- c(300, 500, 750, 1000, 1500),
    g.method <- c("rf", "nnet"),
    gps.method <-
      c("rf&normal" ,
        "linear&boxcox"),
    trimming <- c(-4, 4),
    cov <- c(5),
    method <- c("SR", "CDML", "HI"),
    fold <- c(1, 2, 3, 5),
    responseCurve <- "polynom",
    sd <- 15,
    file <- "polynom.sd15.Rdata" # the file saves data
  )
}
