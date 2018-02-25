


################################################################################
#                                                                              #
#                         Simulation Body Part                                 #
#                                                                              #
################################################################################

#' Simulation Body Part
#'
#' @param kfold number of folds for sample splitting
#' @param data data
#' @param trimLowerBound.pi trimming lower bound of pi estimation
#' @param trimLowerBound.pibar trimming lower bound of pibar estimation
#' @param trimUpperBound.pi trimming upper bound of pi estimation
#' @param trimUpperBound.pibar trimming upper bound of pibar estimation
#' @param method  estimation method we choose from double machine learning method ("CDML")
#'            and simple regression method ("SR"), Hirano & Imbens method ("HI")
#' @param g.method method for regression estimation
#'
#' @param gps.method method for generalized propensity score estimation
#' @param trimLowerBound.t  trimming lower bound of treatment vector
#' @param trimUpperBound.t  trimming upper bound of treatment vector
#' @param model model varies from "IV", "CTE"
#' @param variation  choose from "whole","train.data".
#' Estimate the probability measure by empirical measure
#'                   of covariates. "whole" uses every data sample,
#'                    "train.data" only trainning data for estimating empirical measure.
#'                    As default, we select "whole"
#' @param verbose print in console step by step
#'
#' @return a list of estimated value "yhat" evaluated at "t.grid", with "record" to record details of this simulation
#' @export

simulation <- function(kfold,
                       data = data,
                       trimLowerBound.pi = 0.00001,
                       trimLowerBound.pibar = 0.00001,
                       trimUpperBound.pi = 10,
                       trimUpperBound.pibar = 10,
                       method = "CDML",
                       g.method = "rf",
                       gps.method = "series",
                       trimLowerBound.t = -10,
                       trimUpperBound.t = 10,
                       model,
                       variation = "whole",

                       verbose = TRUE)
{
  ## check whether the number of folds we set is an integer
  check.integer <- function(N) {
    !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
  }

  ## check whether the model name input makes sense
  check.model <- function(name) {
    return(name == "CTE" | name == "IV")
  }

  ## check whether the method name input makes sense
  check.method <- function(name) {
    return(name == "CDML" | name == "SR" | name == "HI")
  }

  ## check whether the variation name is valid
  check.variation <- function(name) {
    return(name == "whole" | name == "data.train")
  }

  ## check whether the g.method name input is valid
  check.g.method <- function(name) {
    return(name == "rf" |
             name == "nnet" |
             name == "trees" | name == "MARS" | name == "lasso")
  }

  ## check whether the gps.method name input is valid
  check.gps.method <- function(name) {
    return(
      name == "series" |
        name == "linear&normal" |
        name == "rf&normal" |
        name == "nnet&normal" | name == "rf&boxcox" |
        name == "nnet&boxcox" |
        name == "linear&boxcox" | name == "boosting&normal" | name == "quantregForest"
    )
  }

  if (!check.model(model)) {
    stop("Please input a valid model name!")
  }

  if (!check.method(method)) {
    stop("Please input a valid method name!")
  }

  if (!check.variation(variation)) {
    stop("Please input a valid variation name!")
  }

  if (!check.integer(kfold)) {
    stop("Please set argument kfold as an integer!")
  }

  if (!check.g.method(g.method)) {
    stop("Please input a valid g.method name!")
  }

  if (!check.gps.method(gps.method)) {
    stop("Please input a valid gps.method name!")
  }

  if (kfold == 1) {
    print("Simulation starts!")

    ## choose which estimation model we choose
    # here we choose double machine learning method
    if (method == "CDML") {
      ## if model is choosen "CTE"
      if (model == "CTE") {
        ## estimation of g and gbar

        object <-
          g.method.estimation(data = data, g.method = g.method)
        g.rf.result <-
          g.predict(
            object = object,
            data.eval = data,
            data.chunk = data,
            mode = method,
            verbose = verbose,
            g.method = g.method
          )

        ## combine g estimation to our data
        ## combine gbar estimation to our data
        gbar <- rowMeans(g.rf.result$matrix)
        data.with.gbar.g <-
          cbind(g = g.rf.result$vector, gbar = gbar, data)


        if (length(gbar) != dim(data)[1] &&
            length(g.rf.result$vector) != dim(data)[1]) {
          stop("g estimation or gbar estimation is invalid!")
        }

        ## estimation of pi, pibar
        object2 <-
          gps.method.estimation(
            data = data,
            gps.method = gps.method,
            treatment.max = trimUpperBound.t,
            treatment.min = trimLowerBound.t,
            verbose = verbose
          )
        pi.series.result <-
          gps.predict(
            data.chunk = data,
            gps.method = gps.method,
            data.eval = data,
            object = object2
          )
        ## combine pi estimation to our data
        ## combine pibar estimation to our data
        pibar <- rowMeans(pi.series.result$matrix)
        data.complete <-
          cbind(pi = pi.series.result$vector,
                pibar = pibar,
                data.with.gbar.g)

        if (length(pibar) != dim(data)[1] &&
            length(pi.series.result$vector) != dim(data)[1]) {
          stop("pi estimation or pibar estimation is invalid!")
        }
      }

      ## if model is choosen "IV"
      if (model == "IV") {
        ## estimation of g and gbar

        ## the instrument does not influence the g estiamtion
        object <-
          g.method.estimation(data = data[, !(names(data) %in% "z")], g.method = g.method)
        g.rf.result <-
          g.predict(
            object = object,
            data.eval = data[, !(names(data) %in% "z")],
            data.chunk = data[, !(names(data) %in% "z")],
            mode = method,
            verbose = verbose,
            g.method = g.method
          )

        ## combine g estimation to our data
        ## combine gbar estimation to our data
        gbar <- rowMeans(g.rf.result$matrix)
        data.with.gbar.g <-
          cbind(g = g.rf.result$vector, gbar = gbar, data)


        if (length(gbar) != dim(data)[1] &&
            length(g.rf.result$vector) != dim(data)[1]) {
          stop("g estimation or gbar estimation is invalid!")
        }

        ## estimation of pi, pibar
        data.iv <- data[, !(names(data) %in% "t")]
        colnames(data.iv)[which(names(data.iv) == "z")] <- "t"

        object2 <-
          gps.method.estimation(
            data = data.iv,
            gps.method = gps.method,
            treatment.min = min(data.iv[, (names(data.iv) %in% "t")]),
            treatment.max = max(data.iv[, (names(data.iv) %in% "t")]),
            verbose = verbose
          )
        pi.series.result <-
          gps.predict(
            data.chunk = data.iv,
            gps.method = gps.method,
            data.eval = data.iv,
            object = object2
          )
        ## combine pi estimation to our data
        ## combine pibar estimation to our data
        pibar <- rowMeans(pi.series.result$matrix)
        data.complete <-
          cbind(pi = pi.series.result$vector,
                pibar = pibar,
                data.with.gbar.g)

        if (length(pibar) != dim(data)[1] &&
            length(pi.series.result$vector) != dim(data)[1]) {
          stop("pi estimation or pibar estimation is invalid!")
        }

      }
    }
    ## choose which estimation model we choose
    # here we choose simple regression method
    if (method == "SR") {

      if(model == "CTE"){
      ## estimation of g an gbar
      object <-
        g.method.estimation(data = data, g.method = g.method)

      ## generate a evaluate grid
      temp.grid <-
        seq(
          from = trimLowerBound.t,
          to = trimUpperBound.t,
          length.out = 100 * (trimUpperBound.t - trimLowerBound.t)
        )

      g.rf.result <-
        g.predict(
          object = object,
          data.eval = temp.grid,
          data.chunk = data,
          mode = method,
          verbose = verbose,
          g.method = g.method

        )

      ## combine gbar estimation to our data
      gbar <- rowMeans(g.rf.result$matrix)
      }

      if(model == "IV"){

        ## estimation of g an gbar
        object <-
          g.method.estimation(data = data[, !(names(data) %in% "z")], g.method = g.method)

        ## generate a evaluate grid
        temp.grid <-
          seq(
            from = trimLowerBound.t,
            to = trimUpperBound.t,
            length.out = 100 * (trimUpperBound.t - trimLowerBound.t)
          )

        g.rf.result <-
          g.predict(
            object = object,
            data.eval = temp.grid,
            data.chunk = data[, !(names(data) %in% "z")],
            mode = method,
            verbose = verbose,
            g.method = g.method

          )

        ## combine gbar estimation to our data
        gbar <- rowMeans(g.rf.result$matrix)
      }

    }

    ## choose which estimation model we choose
    # here we choose Hirano & Imbens method
    if (method == "HI") {
      ## estimation of gps
      object <-
        gps.method.estimation(
          data = data,
          gps.method = gps.method,
          treatment.max = trimUpperBound.t,
          treatment.min = trimLowerBound.t,
          verbose = verbose
        )
      pi.series.result <-
        gps.predict(
          data.chunk = NULL,
          data.eval = data,
          gps.method = gps.method,
          object = object,
          t.grid = seq(
            from = trimLowerBound.t,
            to = trimUpperBound.t,
            length.out = 100 * (trimUpperBound.t - trimLowerBound.t)
          )
        )

      # our REAL complete data to pass down are data.complete and
      # data.complete.comple (complements)
      data.complete.comple <- pi.series.result$matrix
      data.complete <-
        cbind(pi = pi.series.result$vector, data)
    }
  }
  ## we will consider the sample splitting case in below the else bracket
  else{
    print("Simulation starts!")

    ################### data splitting ##########################
    ## data splitting
    split     <- runif(nrow(data))
    split.group   <-
      as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1 / kfold)), include.lowest = TRUE))

    # everything okay, now we store the data in this fold to data.complete
    data.complete <- NULL

    ## focus on individual splitting
    for (jth.fold in 1:kfold) {
      print(paste("fold", jth.fold))

      ii  <- split.group == jth.fold
      nii <- split.group != jth.fold

      # training data & evaluation data
      data.train <- data[nii,]
      data.eval <- data[ii,]

      ## choose which estimation model we choose
      # here we choose double machine learning method
      if (method == "CDML") {
        ## if the chosen model is "CTE"
        if (model == "CTE") {
          ## estimation of g an gbar
          object <-
            g.method.estimation(data = data.train, g.method = g.method)

          if (variation == "data.train") {
            g.rf.result <-
              g.predict(
                object = object,
                data.eval = data.eval,
                data.chunk = data.train,
                verbose = verbose,
                g.method = g.method
              )
          }

          if (variation == "whole") {
            g.rf.result <-
              g.predict(
                object = object,
                data.eval = data.eval,
                data.chunk = data,
                verbose = verbose,
                g.method = g.method
              )
          }
          ## combine g estimation to our data
          ## combine gbar estimation to our data
          gbar <- rowMeans(g.rf.result$matrix)
          data.with.gbar.g <-
            cbind(g = g.rf.result$vector,
                  gbar = gbar,
                  data.eval)

          if (length(gbar) != dim(data.eval)[1] &&
              length(g.rf.result$vector) != dim(data.eval)[1]) {
            stop("g estimation or gbar estimation is invalid!")
          }

          ## estimation of pi, pibar
          object2 <-
            gps.method.estimation(
              data = data.train,
              gps.method = gps.method,
              treatment.max = trimUpperBound.t,
              treatment.min = trimLowerBound.t,
              verbose = verbose
            )
          if (variation == "data.train") {
            pi.series.result <-
              gps.predict(
                data.chunk = data.train,
                data.eval = data.eval,
                gps.method = gps.method,
                object = object2
              )
          }
          if (variation == "whole") {
            pi.series.result <-
              gps.predict(
                data.chunk = data,
                data.eval = data.eval,
                gps.method = gps.method,
                object = object2
              )
          }
          ## combine pi estimation to our data
          ## combine pibar estimation to our data
          pibar <- rowMeans(pi.series.result$matrix)

          data.complete.foldwise <-
            cbind(pi = pi.series.result$vector,
                  pibar = pibar,
                  data.with.gbar.g)

          if (length(pibar) != dim(data.eval)[1] &&
              length(pi.series.result$vector) != dim(data.eval)[1]) {
            stop("pi estimation or pibar estimation is invalid!")
          }


          data.complete <-
            rbind(data.complete, data.complete.foldwise)
        }

        ## if the chosen model is "IV"
        if (model == "IV") {
          ## estimation of g an gbar
          object <-
            g.method.estimation(data = data.train[, !(names(data.train) %in% "z")], g.method = g.method)

          if (variation == "data.train") {
            g.rf.result <-
              g.predict(
                object = object,
                data.eval = data.eval[, !(names(data.train) %in% "z")],
                data.chunk = data.train[, !(names(data.train) %in% "z")],
                verbose = verbose,
                g.method = g.method
              )
          }

          if (variation == "whole") {
            g.rf.result <-
              g.predict(
                object = object,
                data.eval = data.eval[, !(names(data.train) %in% "z")],
                data.chunk = data[, !(names(data.train) %in% "z")],
                verbose = verbose,
                g.method = g.method
              )
          }
          ## combine g estimation to our data
          ## combine gbar estimation to our data
          gbar <- rowMeans(g.rf.result$matrix)
          data.with.gbar.g <-
            cbind(g = g.rf.result$vector,
                  gbar = gbar,
                  data.eval)

          if (length(gbar) != dim(data.eval)[1] &&
              length(g.rf.result$vector) != dim(data.eval)[1]) {
            stop("g estimation or gbar estimation is invalid!")
          }

          ## estimation of pi, pibar

          data.iv <- data[, !(names(data) %in% "t")]
          colnames(data.iv)[which(names(data.iv) == "z")] <- "t"

          data.train.iv <-
            data.train[, !(names(data.train) %in% "t")]
          colnames(data.train.iv)[which(names(data.train.iv) == "z")] <-
            "t"

          data.eval.iv <- data.eval[, !(names(data.eval) %in% "t")]
          colnames(data.eval.iv)[which(names(data.eval.iv) == "z")] <-
            "t"

          object2 <-
            gps.method.estimation(
              data = data.train.iv,
              gps.method = gps.method,
              treatment.min = min(data.iv[, (names(data.iv) %in% "t")]),
              treatment.max = max(data.iv[, (names(data.iv) %in% "t")]),
              verbose = verbose
            )
          if (variation == "data.train") {
            pi.series.result <-
              gps.predict(
                data.chunk = data.train.iv,
                data.eval = data.eval.iv,
                gps.method = gps.method,
                object = object2
              )
          }
          if (variation == "whole") {
            pi.series.result <-
              gps.predict(
                data.chunk = data.iv,
                data.eval = data.eval.iv,
                gps.method = gps.method,
                object = object2
              )
          }
          ## combine pi estimation to our data
          ## combine pibar estimation to our data
          pibar <- rowMeans(pi.series.result$matrix)

          data.complete.foldwise <-
            cbind(pi = pi.series.result$vector,
                  pibar = pibar,
                  data.with.gbar.g)
          if (length(pibar) != dim(data.eval)[1] &&
              length(pi.series.result$vector) != dim(data.eval)[1]) {
            stop("pi estimation or pibar estimation is invalid!")
          }


          data.complete <-
            rbind(data.complete, data.complete.foldwise)

        }
      }



      #################### the end of data sampling part ####################
      cat("\n")
    }
  }

  #####################       post-processing ######################

  ## post-processing for CDML method
  if (method == "CDML") {
    ####### we have to filter the conditional density so that ub > pi, pibar > lb ######
    data.complete <-
      data.complete[(
        data.complete$pibar > trimLowerBound.pibar &
          trimUpperBound.pibar > data.complete$pibar
      ),]
    data.complete <-
      data.complete[(data.complete$pi > trimLowerBound.pi &
                       trimUpperBound.pi > data.complete$pi),]

    row.names(data.complete) <- 1:nrow(data.complete)
    ##############################################################

    ## The formula: we use here CDML
    pseudo.out <-
      (data.complete$y - data.complete$g) / (data.complete$pi / data.complete$pibar) + data.complete$gbar
    data.complete <- cbind(data.complete, ybar = pseudo.out)

    ## to what degree out debiased technique will influence the pseudo-output
    influence.pseudo.out <-
      abs(
        (data.complete$y - data.complete$g) / (data.complete$pi / data.complete$pibar) /
          data.complete$gbar
      )

    ###################### local linear regression part! #########################
    treatment <- data.complete$t
    response <- data.complete$ybar

    ## centralize and plot!
    x <- quantile(response, c(0.01, 0.99))
    response.cleaned <-
      response[response >= x[1] & response <= x[2]]
    treatment.cleaned <-
      treatment[response >= x[1] & response <= x[2]]
    #plot(treatment.cleaned, response.cleaned)

    h <- 0
    # choose bandwidth by https://stat.ethz.ch/R-manual/R-devel/library/KernSmooth/html/dpill.html
    try(h <- dpill(treatment, response))
    # this happens sometimes that by dpill we cannot select bandwidth h
    # hence below try to fix this scenerio and report when h is NA

    if (is.na(h) | h <= 0) {
      print("Caution: dpill fails to choose h. We will choose h manually")

      # manually select bandwidth from a given range
      bw.from = 0.01 # bandwidth has to be strictly positive
      bw.to = (max(treatment) - min(treatment)) / 3 # bandwidth cant be too wide
      bw.length.out = 100
      epsGrid = seq(bw.from, bw.to, length.out = bw.length.out) # range of tuning parameters from the kernel (bandwidth)

      error = rep(NA, length(epsGrid)) # vector to store the errors of the models
      # when using different parameter value

      # tuning epsilon:
      for (ii in 1:length(epsGrid))
      {
        eps = epsGrid[ii]
        object = locpoly(treatment,
                         response,
                         bandwidth = eps,
                         degree = 1)
        error[ii] = sum((object$y - response) ^ 2)
      }

      bestEps = epsGrid[which.min(error)] # best bandwidth

      # run this seriesCDE again with best bandwidth and set to choose Delta
      final.result =  locpoly(
        treatment,
        response,
        bandwidth = bestEps,
        degree = 1,
        range.x = c(trimLowerBound.t, trimUpperBound.t),
        gridsize = 100 * (trimUpperBound.t - trimLowerBound.t)
      )
      #lines(final.result)
      print(h)
    } else{
      # final.result as the result of local linear regression
      final.result <-
        locpoly(
          treatment,
          response,
          bandwidth = h,
          degree = 1,
          range.x = c(trimLowerBound.t, trimUpperBound.t),
          gridsize = 100 * (trimUpperBound.t - trimLowerBound.t)
        )
      #lines(final.result)
    }

    ########## NA-value of CDML estimation filling #############
    if (any(is.na(final.result$y)) |
        any(is.infinite(final.result$y))) {
      ## to inform the user that we have done this NA-value filling process
      print("CDML estimation cannot estimate some values at the edge.")
      cat("\n")
      print("Hence we fill the missing edge values with the closest value to the missing ones.")
      print("Or we change the Inf value to the closest non-Inf value")

      # this problem is because of locpoly does weirdly when estimating treatment effects far from data
      ### convert value in result.cdml$yhat to NA if the value == infinite
      indices <- which(is.infinite(final.result$y))

      ##debug
      print(indices)

      for (indice in indices) {
        final.result$y[indice] <- NA
      }
      ### change the value in result.cdml$yhat to the closest value, if the value == NA
      indices <- which(is.na(final.result$y))

      ## debug
      print(indices)
      for (indice in indices) {
        if (abs(indice - 1) > abs(indice - length(final.result$y))) {
          final.result$y[indice] <-
            final.result$y[max(which(!is.na(final.result$y)))]
        } else{
          final.result$y[indice] <-
            final.result$y[min(which(!is.na(final.result$y)))]
        }
      }
    }

  }

  ## post-processing for simple regression method
  # if (method == "SR") {
  # }

  ## post-processing for Hirano & Imbens method
  if (method == "HI") {
    hi.formula <- y ~ pi + I(pi ^ 2) + t + I(t ^ 2) + t * pi
    hi.model <- lm(hi.formula, data = data.complete)


    t.grid <-
      seq(
        from = trimLowerBound.t,
        to = trimUpperBound.t,
        length.out = 100 * (trimUpperBound.t - trimLowerBound.t)
      )
    yhat <-
      rep(NA, length(t.grid)) # vector to save the final result

    for (i in 1:length(t.grid)) {
      data.temp <- data.complete
      data.temp$t <- t.grid[i]
      data.temp$pi <- data.complete.comple[i, ]
      ## do things!
      yhat[i] <- mean(predict(hi.model, newdata = data.temp))
    }

    #plot(x = t.grid,
        # y = yhat,
        # type = "l")
  }

  ################# returns ###################
  ########### legend to record ##########

  ## record the number of covariates
  if (model == "CTE") {
    numofcov <- dim(data)[2] - 2
  }

  if (model == "IV") {
    numofcov <- dim(data)[2] - 3
  }


  legend <-
    list(
      "g.method" = g.method,
      "gps.method" = gps.method,
      "trimming" = c(trimLowerBound.t, trimUpperBound.t),
      "method" = method,
      "fold" = kfold,

      "covariates" = numofcov,
      "samples" = dim(data)[1]
    )
  # ######################################
  if (method == "HI") {
    return.list <-
      list("yhat" = yhat,
        #   "t.grid" = t.grid,
           "record" = legend)
  }

  if (method == "CDML") {
    return.list <-
      list("yhat" = final.result$y,
         #  "t.grid" = final.result$x,
           "record" = legend)
  }

  if (method == "SR") {
    return.list <-
      list("yhat" = gbar,
        #   "t.grid" = temp.grid,
           "record" = legend)
    #plot(return.list$t.grid, return.list$yhat, type = "l")
  }
  return(return.list)
}
