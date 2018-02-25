#' The alternative of running simulation
#'
#' @return file saved with name specified
#' @export
#'
simulation.alternative <- function() {
  ncores <-
    min(detectCores(all.tests = FALSE, logical = TRUE), 44)

  file <- "simu.alter.Rdata"

  rep <- 7
  simu.times <- 2
  model <- "CTE"
  trimming <- c(-4, 4)


  method <- c("CDML")
  g.method <- c("rf", "nnet")
  gps.method <- c("rf&normal")
  fold <- c(1, 2)

  #################### data.generation set up ###################


  sample.size.option <- c(200, 300)
  num.cov.option <- c(2, 3, 5, 10)
  treatment.sd.option <- c(1, 2, 3, 5, 8, 10)
  struct.model.option <-
    c("linear", "polynom", "polynom2", "polynom3")


  sample.size.generator <-
    sample(sample.size.option, simu.times, replace = TRUE)
  num.cov.generator <-
    sample(num.cov.option, simu.times, replace = TRUE)
  treatment.sd.generator <-
    sample(treatment.sd.option, simu.times, replace = TRUE)
  struct.model.generator <-
    sample(struct.model.option, simu.times, replace = TRUE)

  ####################### simulation.set.up ##################################

  simu.setups.inbetween1 <- NULL
  simu.setups.inbetween2 <- NULL
  simu.setups.inbetween3 <- NULL

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

  simu.setups <- simu.setups.inbetween

  ######################## indicator of time to wait ##################################

  print("fire, taget number:")
  cat("\n")
  print(simu.times * length(simu.setups$method))

  ####################### data.list.list ##################################

  data.list.list <-
    mapply(
      simu = rep,
      samples = sample.size.generator,
      covariates = num.cov.generator,
      parametricCurveOption = struct.model.generator,
      sd = treatment.sd.generator,
      FUN = data.generation.wrap,
      MoreArgs = list(
        trimLowerBound = trimming[1],
        trimUpperBound = trimming[2],
        model = model,
        ncores = ncores
      ),
      SIMPLIFY = FALSE
    )

  ################## true.data.list ###################################

  true.data.list <- mcmapply(
    FUN = true.data.generation,
    parametricCurveOption <- struct.model.generator,
    MoreArgs = list(
      numOfCovariates = 2,
      trimLowerBound = trimming[1],
      trimUpperBound = trimming[2],
      distributionForCovariates = "normal",
      distributionForNoise = "normal",
      distributionForTreatment = "normal",
      nrep = 1000000
    ),
    SIMPLIFY = FALSE,
    mc.cores = ncores
  )

  ################# result.list ###################

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

  ################# save~1 ####################
  save(data.list.list,
       simu.setups,
       true.data.list,
       res.list.list,
       file = file)

  ######################## error #############################

  error.list.list <- mapply(
    result.list = res.list.list,
    true.data = true.data.list,
    FUN = err.w,
    MoreArgs = list(ncores = ncores),
    SIMPLIFY = FALSE
  )

  bias.list <-
    bias.evaluation.wrap(res.list.list = res.list.list, true.data = true.data.list)

  ##################### save~2 ##############################

  save(
    data.list.list,
    simu.setups,
    true.data.list,
    res.list.list,
    error.list.list,
    bias.list,
    file = file
  )
}
