






################################################################################
#                                                                              #
#                            Demostration                                      #
#                                                                              #
################################################################################



demo <- function(model = "CTE",
                 fold = 5,
                 g.method,
                 gps.method,
                 numOfSamples = 200) {
  #set.seed(101) # for reproducibility

  curve <- "polynom"

  data <- data.generation(
    model = model,
    numOfSamples = numOfSamples,
    numOfPreSamples = NULL,
    numOfCovariates = 5,
    distributionForCovariates = "normal",
    distributionForIVs = "normal",
    distributionForNoise = "normal",
    distributionForTreatment = "normal",
    parametricCurveOption = curve,
    trimLowerBound = -4,
    trimUpperBound = 4,
    sd = 8
  )

  true.data <-
    true.data.generation(
      parametricCurveOption = curve,
      trimLowerBound = -4,
      trimUpperBound = 4
    )

  print("Data Generating ...")
  Sys.sleep(2)

  plot(data$t,
       data$y,
       xlab = "treatment",
       ylab = "response")

  print("Underlying true curve generating ...")
  Sys.sleep(2)

  lines(true.data$grid,
        true.data$true.y,
        col = "red",
        lwd = 3)


  result <- simulation(
    kfold = fold ,
    data = data,
    method = "CDML",
    g.method = g.method,
    gps.method = gps.method,
    trimLowerBound.t = -4,
    trimUpperBound.t = 4,
    model = model,
    verbose = FALSE
  )

  lines(true.data$grid,
        result$yhat,
        col = "blue",
        lwd = 3)

  result2 <- simulation(
    kfold = 1 ,
    data = data,
    method = "SR",
    g.method = g.method,
    gps.method = gps.method,
    trimLowerBound.t = -4,
    trimUpperBound.t = 4,
    model = model,
    verbose = FALSE
  )

  lines(true.data$grid,
        result2$yhat,
        col = "cyan",
        lwd = 3)

  result3 <- simulation(
    kfold = 1 ,
    data = data,
    method = "HI",
    g.method = g.method,
    gps.method = gps.method,
    trimLowerBound.t = -4,
    trimUpperBound.t = 4,
    model = model,
    verbose = FALSE
  )

  lines(true.data$grid,
        result3$yhat,
        col = "purple",
        lwd = 3)

}
