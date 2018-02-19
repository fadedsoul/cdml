################################################################################
#                            model: IV/CTE model                               #
#           response.y.: polynom 1,2,3 and linear case                         #
#          g.method <- c("nnet", "rf","trees", "MARS")                         #
#     gps.method <- c( "series", "boosting&normal", "linear&boxcox")           #
################################################################################

################ FUNCTION ###########################
application.samples <- function(model) {

  if (model == "IV") {
    file1 <- "linear.sample.iv.Rdata"
    file2 <- "polynom.sample.iv.Rdata"
    file3 <- "polynom2.sample.iv.Rdata"
    file4 <- "polynom3.sample.iv.Rdata"
  }
  if(model == "CTE") {
    file1 <- "linear.sample.cte.Rdata"
    file2 <- "polynom.sample.cte.Rdata"
    file3 <- "polynom2.sample.cte.Rdata"
    file4 <- "polynom3.sample.cte.Rdata"
  }

 running.simulation(
   ########### enter all parameters needed for running simulations ##########
   model <- model,
   simu <- 40,
   samples <- c(200, 300, 500, 750, 1000, 1500, 2000),
   g.method <- c("rf","nnet","trees"),
   gps.method <-
     c("series",
      "rf&normal" ,
      "linear&boxcox"
      ),
   trimming <- c(-4, 4),
   cov <- c(5),
   method <- c("SR","CDML","HI"),
   fold <- c(1,2,3,5),
   responseCurve <- "linear",
   file <- file1 # the file saves data
 )

running.simulation(
  ########### enter all parameters needed for running simulations ##########
  model <- model,
  simu <- 40,
  samples <- c(200, 300, 500, 750, 1000, 1500, 2000),
  g.method <- c("rf","nnet","trees"),
  gps.method <-
    c("series",
      "rf&normal" ,
      "linear&boxcox"
    ),
  trimming <- c(-4, 4),
  cov <- c(5),
  method <- c("SR","CDML","HI"),
  fold <- c(1,2,3,5),
  responseCurve <- "polynom",
  file <- file2 # the file saves data
)

running.simulation(
  ########### enter all parameters needed for running simulations ##########
  model <- model,
  simu <- 40,
  samples <- c(200, 300, 500, 750, 1000, 1500, 2000),
  g.method <- c("rf","nnet","trees"),
  gps.method <-
    c("series",
      "rf&normal" ,
      "linear&boxcox"
    ),
  trimming <- c(-4, 4),
  cov <- c(5),
  method <- c("SR","CDML","HI"),
  fold <- c(1,2,3,5),
  responseCurve <- "polynom2",
  file <- file3 # the file saves data
)

running.simulation(
  ########### enter all parameters needed for running simulations ##########
  model <- model,
  simu <- 40,
  samples <- c(200, 300, 500, 750, 1000, 1500, 2000),
  g.method <- c("rf","nnet","trees"),
  gps.method <-
    c("series",
      "rf&normal" ,
      "linear&boxcox"
    ),
  trimming <- c(-4, 4),
  cov <- c(5),
  method <- c("SR","CDML","HI"),
  fold <- c(1,2,3,5),
  responseCurve <- "polynom3",
  file <- file4 # the file saves data
)
}


##########################################################
#                                                        #
#    Remark: the plotting function is saved locally !    #
#          Ask if needed.                                #
##########################################################
