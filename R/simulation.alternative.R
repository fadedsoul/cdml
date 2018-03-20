#' The alternative of running simulation
#'
#' @param model model varies from "IV", "CTE"
#' @param rep within each set-up, number of repetition that will be done
#' @param file file that we save our return into
#' @param method estimation method we choose from double machine learning method ("CDML")
#'            and simple regression method ("SR"), Hirano & Imbens method ("HI")
#' @param simu.times number of simulations one would run
#' @param g.method a vector of method for regression estimation
#' @param gps.method a vector of method for generalized propensity score estimation
#' @param fold number of folds for sample splitting
#' @param sample.size.option a vector of number of samples that will be selected from
#' @param num.cov.option a vector of number of covariates that will be selected from
#' @param treatment.sd.option a vector of standard deviation of treatment generation that will be selected from
#' @param struct.model.option a vector of parametric curves of response generation, typically choose from "linear", "polynom", "polynom2", "polynom3", "mixture".
#' @param detoured if FALSE method will be done normally;
#'             otherwise a special simulation will be examined, where CDML is chosen method, CTE is chosen model,
#'             gps.method "series" is not allowed.
#' @param keep.same in the case detoured is true, whether we use the same gps.method to estimate T based on X
#' @return file saved with name specified
#' @export
#'
simulation.alternative <-
  function(model = "CTE",
           rep = 7,
           file =
             #"simu.alter.Rdata",
             "simu.alter.rds",
           method = c("CDML"),
           simu.times = 2,
           g.method = c("rf", "nnet"),
           gps.method = c("linear&normal"),
           fold = c(1),
           sample.size.option = c(200, 300,500, 800, 1000, 1500,2000),
           num.cov.option = c(1,2, 3, 5, 8,10,15, 20),
           treatment.sd.option = c(0.5, 1, 2, 3, 5, 8, 10),
           struct.model.option =
             c("linear", "polynom", "polynom2", "polynom3"),
           detoured = FALSE,
           keep.same = TRUE) {
    ###################### Pre-Set-Up ######################

    set.seed(100)

    ncores <-
      min(detectCores(all.tests = FALSE, logical = TRUE), 44)

    trimming <- c(-4, 4)

    if (model == "IV") {
      if (detoured == TRUE) {
        error(
          "detoured has to be set false if you are using this function for Instrumental Variable model!"
        )
      }
    }

    #################### data.generation set up ###################

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
      ) == FALSE), ]
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

    ############ Normally Speaking ##################
    if (detoured == FALSE) {

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
      ### just to be safe
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
    }

    ################# if we are about to use the detoured trick ################
    if(detoured == TRUE){

      if(keep.same == TRUE){

        dat.list.list <- mapply(
          data = data.list.list,
          FUN = simulation.wrap,
          MoreArgs = list(
            trimLowerBound.t = trimming[1],
            trimUpperBound.t = trimming[2],
            model = model,
            ncores = ncores,
            detoured = TRUE,
            method = "CDML",
            kfold = fold,
            g.method = "rf",
            gps.method = gps.method
          ),
          SIMPLIFY = FALSE
        )
      }

      res.list.list <- mapply(
        data = rep(dat.list.list, each = length(simu.setups$method)),
        method = rep(simu.setups$method, times = length(data.list.list)),
        kfold = rep(simu.setups$fold, times = length(data.list.list)),
        g.method = rep(simu.setups$g.method, times = length(data.list.list)),
        gps.method = rep(simu.setups$gps.method, times = length(data.list.list)),
        FUN = simulation.wrap,
        MoreArgs = list(
          trimLowerBound.t = trimming[1],
          trimUpperBound.t = trimming[2],
          model = "IV",
          ncores = ncores
        ),
        SIMPLIFY = FALSE
      )

   saveRDS(list(data.list.list = dat.list.list,
                res.list.list = res.list.list), file = file)

   error.list.list <- mapply(
     result.list = res.list.list,
     true.data = true.data.list,
     FUN = err.w,
     MoreArgs = list(ncores = ncores),
     SIMPLIFY = FALSE
   )

   bias.list <-
     bias.evaluation.wrap(res.list.list = res.list.list, true.data = true.data.list)

   saveRDS(
     list(
       data.list.list = data.list.list,
       simu.setups = simu.setups,
       true.data.list = true.data.list,
       res.list.list = res.list.list,
       error.list.list = error.list.list,
       bias.list = bias.list
     ),
     file = file
   ) ## << convention use file "<something>.rds"

    }
    ##################### save~2 ##############################

    if(detoured == FALSE){
    ## Alternative  to  save()  and  load()
    saveRDS(
      list(
        data.list.list = data.list.list,
        simu.setups = simu.setups,
        true.data.list = true.data.list,
        res.list.list = res.list.list,
        error.list.list = error.list.list,
        bias.list = bias.list
      ),
      file = file
    ) ## << convention use file "<something>.rds"
    }
  }
