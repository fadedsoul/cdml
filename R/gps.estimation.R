

################################################################################
#                                                                              #
#                         gps.method.estimation                                #
#                                                                              #
################################################################################

#' Generalized Propensity Score Estimation
#'
#' @param data data
#' @param gps.method methods for generalized propensity score estimation, choose from "series", "linear&normal", "boxcox"
#' @param n.prop default = 0.7 meaning 70\% for training usage, 30\% for validation usage in method "series"
#' @param bw.length.out number of tested points of tuning parameters from the kernel (bandwidth) in method "series"
#' @param bw.from lower bound of tuning parameters from the kernel (bandwidth) in method "series"
#' @param bw.to upper bound of tuning parameters from the kernel (bandwidth) in method "series"
#' @param treatment.min  in what interval conditional density shall be evaluated
#' @param treatment.max  in what interval conditional density shall be evaluated
#' @param verbose print in console the procedure step by step
#'
#' @return object of generalized propensity score estimation
gps.method.estimation <- function(data,
                                  gps.method = "series",
                                  n.prop = 0.7,
                                  bw.length.out = 10,
                                  bw.from = 0.5,
                                  bw.to = 9.5,
                                  treatment.min = -10,
                                  treatment.max = 10,
                                  verbose = TRUE)

{
  if (!("t" %in% names(data) && "y" %in% names(data))) {
    stop("The input data is invalid!")
  } # check whether the data is written in a legal way

  #### SeriesSpecCDE ####
  if (gps.method == "series") {
    n <- dim(data)[1] # number of samples
    indices <-
      sample(n, n * n.prop, replace = FALSE)

    data.train <-
      data[indices,][,!(names(data) %in% c("y"))] # training data
    data.valid <-
      data[-indices,][,!(names(data) %in% c("y"))] # validation data

    nXMax <-
      100 # Maximum possible number of components
    # of the series expansion in the X direction

    repeat {
      # tune parameters associated to kernel, and I (nZBest) and J (nXBest)
      epsGrid = seq(bw.from, bw.to, length.out = bw.length.out) # range of tuning parameters from the kernel (bandwidth)

      error = rep(NA, length(epsGrid)) # vector to store the errors of the models
      # when using different parameter value

      # tuning epsilon:
      for (ii in 1:length(epsGrid))
      {
        #print(paste(ii / length(epsGrid) * 100, "percent complete", sep = ""))
        eps = epsGrid[ii]
        object = specSeriesCDE::seriesCDE(
          xTrain = data.matrix(data.train[,!(names(data.train) %in% c("t"))]),
          xValidation = data.matrix(data.valid[,!(names(data.valid) %in% c("t"))]),
          zTrain = data.train[, (names(data.train) %in% c("t"))],
          zValidation = data.valid[, (names(data.valid) %in% c("t"))],
          kernelFunction = specSeriesCDE::radialKernel,
          extraKernel = list("eps.val" = eps),
          chooseDelta = FALSE,
          nXMax = nXMax,
          verbose = verbose
        )
        error[ii] = object$bestError
        #cat("\n")
      }

      bestEps = epsGrid[which.min(error)] # best bandwidth

      # run this seriesCDE again with best bandwidth and set to choose Delta
      object = specSeriesCDE::seriesCDE(
        xTrain = data.matrix(data.train[,!(names(data.train) %in% c("t"))]),
        xValidation = data.matrix(data.valid[,!(names(data.valid) %in% c("t"))]),
        zTrain = data.train[, (names(data.train) %in% c("t"))],
        zValidation = data.valid[, (names(data.valid) %in% c("t"))],
        kernelFunction = specSeriesCDE::radialKernel,
        extraKernel = list("eps.val" = bestEps),
        chooseDelta = TRUE,
        zMin = treatment.min,
        zMax = treatment.max,
        nXMax = nXMax,
        verbose = verbose
      )

      # display chosen parameters
      print(object$nXBest)

      if (object$nXBest < 0.9 * nXMax) {
        break
      } else{
        print("nXBest is close to nXMax, we now increase nXMax")

        nXMax = 4 * object$nXBest # we choose a larger nXMax and run the whole process again
      }
    }

  }

  #### linear model of t,x, and the noise is normal distributed ####
  if (gps.method == "linear&normal") {
    # filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    lm.res <- lm(t ~ ., data = data)

    mean <- mean(lm.res$residuals)
    sd <- sd(lm.res$residuals)

    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(lm.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
  }

  #### rf model of t,x, and the noise is normal distributed ####
  if (gps.method == "rf&normal") {
    # filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    rf.res <- randomForest(t ~ ., data = data, ntrees = 1000)

    mean <- mean(data$t - rf.res$predicted)
    sd <- sd(data$t - rf.res$predicted)

    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(rf.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
  }

  #### nnet model of t,x, and the noise is normal distributed ####
  if (gps.method == "nnet&normal") {
    # filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    nnet.res <- nnet(t ~ .,
                     size = 8,
                     linout = TRUE,
                     data = data)

    mean <- mean(nnet.res$residuals)
    sd <- sd(nnet.res$residuals)

    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(nnet.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
  }

  #### boosting model of t,x and the noise is normal distributed ####
  if (gps.method == "boosting&normal") {
    # filter out the y column
    data <- data[,!(names(data) %in% c("y"))]

    ### A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments.
    ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4749263/pdf/nihms725897.pdf
    F.aac.iter = function(i,
                          data,
                          ps.model,
                          ps.num,
                          rep,
                          criterion) {
      # i: number of iterations (trees)
      # data: dataset containing the treatment and the covariates
      # ps.model: the boosting model to estimate p(T_i|X_i)
      # ps.num: the estimated p(T_i)
      # rep: number of replications in bootstrap
      # criterion: the correlation metric used as the stopping criterion
      GBM.fitted = predict(
        ps.model,
        newdata = data,
        n.trees = floor(i),
        type = "response"
      )
      ps.den = dnorm((data$t - GBM.fitted) / sd(data$t - GBM.fitted), 0, 1)
      wt = ps.num / ps.den
      aac_iter = rep(NA, rep)
      for (i in 1:rep) {
        bo = sample(1:dim(data)[1], replace = TRUE, prob = wt)
        newsample = data[bo, ]
        j.drop = match(c("t"), names(data))
        j.drop = j.drop[!is.na(j.drop)]
        x = newsample[, -j.drop]
        if (criterion == "spearman" | criterion == "kendall") {
          ac = apply(
            x,
            MARGIN = 2,
            FUN = cor,
            y = newsample$t,
            method = criterion
          )
        } else if (criterion == "distance") {
          ac = apply(x,
                     MARGIN = 2,
                     FUN = dcor,
                     y = newsample$t)
        } else if (criterion == "pearson") {
          ac = matrix(NA, dim(x)[2], 1)
          for (j in 1:dim(x)[2]) {
            ac[j] = ifelse (
              !is.factor(x[, j]),
              cor(newsample$t, x[, j],
                  method = criterion),
              polyserial(newsample$t, x[, j])
            )
          }
        } else
          print("The criterion is not correctly specified")
        aac_iter[i] = mean(abs(1 / 2 * log((1 + ac) / (1 - ac))), na.rm =
                             TRUE)
      }
      aac = mean(aac_iter)
      return(aac)
    }
    # Find the optimal number of trees using Pearson/polyserial correlation
    if (shapiro.test(data$t)$p.value < 0.1) {
      model.num = lm(t ~ 1, data = data)
      ps.num = dnorm((data$t - model.num$fitted) / (summary(model.num))$sigma, 0, 1)
    }
    else{
      model.num <- npudens(data$t)
      ps.num <- predict(model.num, newdata = data$t)
    }
    model.den = gbm(
      t ~ .,
      data = data,
      shrinkage = 0.0005,
      interaction.depth = 4,
      distribution = "gaussian",
      n.trees = 20000
    )
    opt = optimize(
      F.aac.iter,
      interval = c(1, 20000),
      data = data,
      ps.model = model.den,
      ps.num = ps.num,
      rep = 50,
      criterion = "pearson"
    )
    best.aac.iter = opt$minimum
    best.aac = opt$objective

    res <-  data$t - predict(model.den, newdata = data,   n.trees = floor(best.aac.iter),
                                   type = "response")
    sd <- sd(res)
    mean <- mean(res)

    ## returns this object as usual as a function
    object <- function(data.eval) {
      fitted = predict(
        model.den,
        newdata = data.eval,
        n.trees = floor(best.aac.iter),
        type = "response"
      )
      return(dnorm((data.eval$t - fitted) ,
                 mean = mean, sd = sd) )
    }
  }

  ####neural network box-cox ####
  if (gps.method == "nnet&boxcox") {
    ## filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    nnet.res <- nnet(t ~ .,
                     size = 8,
                     linout = TRUE,
                     data = data)
    nnet.res.res <- nnet.res$residuals
    nnet.res.exp <- exp(nnet.res.res)

    lambda <-
      EnvStats::boxcox(as.vector(nnet.res.exp),
                       optimize = TRUE,
                       eps = 0.01)$lambda

    if (lambda != 0) {
      mean <- mean((nnet.res.exp ^ lambda - 1) / lambda)
      sd <- sd((nnet.res.exp ^ lambda - 1) / lambda)
    } else {
      mean <- mean(nnet.res.res)
      sd <- sd(nnet.res.res)
    }

    if (lambda != 0) {
      object <- function(data.eval) {
        return(dnorm((exp(
          data.eval$t - predict(nnet.res, data.eval)
        ) ^ lambda - 1) / lambda,
        mean = mean,
        sd = sd,
        log = FALSE
        ))
      }
    } else{
      object <- function(data.eval) {
        return(dnorm(
          data.eval$t - predict(nnet.res, data.eval),
          mean = mean,
          sd = sd,
          log = FALSE
        ))
      }
    }
  }

  #### random forest box-cox ####
  if (gps.method == "rf&boxcox") {
    ## filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    rf.res <- randomForest(t ~ ., data = data, ntrees = 500)
    rf.res.res <- data$t - rf.res$predicted
    rf.res.exp <- exp(rf.res.res)

    lambda <-
      EnvStats::boxcox(as.vector(rf.res.exp),
                       optimize = TRUE,
                       eps = 0.01)$lambda

    if (lambda != 0) {
      mean <- mean((rf.res.exp ^ lambda - 1) / lambda)
      sd <- sd((rf.res.exp ^ lambda - 1) / lambda)
    } else {
      mean <- mean(rf.res.res)
      sd <- sd(rf.res.res)
    }

    if (lambda != 0) {
      object <- function(data.eval) {
        return(dnorm((exp(
          data.eval$t - predict(rf.res, data.eval)
        ) ^ lambda - 1) / lambda,
        mean = mean,
        sd = sd,
        log = FALSE
        ))
      }
    } else{
      object <- function(data.eval) {
        return(dnorm(
          data.eval$t - predict(rf.res, data.eval),
          mean = mean,
          sd = sd,
          log = FALSE
        ))
      }
    }
  }

  #### box-cox  linear ####
  if (gps.method == "linear&boxcox") {
    ## filter out the y column
    data <- data[,!(names(data) %in% c("y"))]
    lm.res <- lm(t ~ ., data = data)
    lm.res.res <- lm.res$residuals
    lm.res.exp <- exp(lm.res.res)

    lambda <-
      EnvStats::boxcox(as.vector(lm.res.exp),
                       optimize = TRUE,
                       eps = 0.01)$lambda

    if (lambda != 0) {
      mean <- mean((lm.res.exp ^ lambda - 1) / lambda)
      sd <- sd((lm.res.exp ^ lambda - 1) / lambda)
    } else {
      mean <- mean(lm.res.res)
      sd <- sd(lm.res.res)
    }

    if (lambda != 0) {
      object <- function(data.eval) {
        return(dnorm((exp(
          data.eval$t - predict(lm.res, data.eval)
        ) ^ lambda - 1) / lambda,
        mean = mean,
        sd = sd,
        log = FALSE
        ))
      }
    } else{
      object <- function(data.eval) {
        return(dnorm(
          data.eval$t - predict(lm.res, data.eval),
          mean = mean,
          sd = sd,
          log = FALSE
        ))
      }
    }
  }

  return(object)
}

################################################################################
#                                                                              #
#                                gps.predict                                   #
#                                                                              #
################################################################################

#' Generalize Propensity Score Prediction
#'
#' @param object an object from gps.method.estimation function
#' @param gps.method selected generalized propensity score method
#' @param data.chunk whole data, covariates are of interest. defualt = NULL
#' @param data.eval data to be evaluated, the treatment is of interest
#' @param t.grid alternate parameter to set up if data.chunk is NULL
#' @param grid.length how dense the evaluation point should be like
#' @param verbose print in console the procedure step by step
#'
#' @return a vector and a matrix
#' For the matrix elements, they stand for the predicted values of:
#' row means given certain treatment from data.eval
#' column means given certain covariates from data.chunk.
#' For the vector: we list in order the predicted value as \eqn{i=1,...,n} for \eqn{(T_i, X_i)}
#'
#'
gps.predict <-
  function(object,
           gps.method = "series",
           data.chunk = NULL,
           data.eval,
           t.grid = NULL,
           grid.length = 1000,
           verbose = FALSE)
  {
    # the value to be returned later
    result.mat <- NULL
    result.vec <- NULL

    # if data.chunk variable is degenerated, then we skip the matrix-result calculation
    if (!is.null(data.chunk)) {
      result.mat <-
        matrix(nrow = dim(data.eval)[1],
               ncol = dim(data.chunk)[1]) # matrix to save the final result
      # row means given certain treatment from data.eval
      # column means given certain covariates from data.complete

      ### series method ###
      if (gps.method == "series") {
        ## return matrix saving pi(T|x) for various x, used in \int_\CX \pi(T|x)dP_x
        x.data.chunk <-
          data.matrix(data.chunk[,!(names(data.chunk) %in% c("t", "y"))])
        result.raw <-
          specSeriesCDE::predictCDE(object = object,
                     xTest = x.data.chunk,
                     B = grid.length)

        ## save the conditional density value into the result.mat
        for (ii in 1:dim(data.eval)[1]) {
          treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]

          ind.1 <- which.min(abs(result.raw$z - treatmentValue))

          if (ind.1 == 1) {
            ind.2 <- ind.1 + 1
          } else{
            if (ind.1 == grid.length) {
              ind.2 <- ind.1 - 1
            } else{
              if (abs(result.raw$z[ind.1 + 1] - treatmentValue) <
                  abs(result.raw$z[ind.1 - 1] - treatmentValue)) {
                ind.2 <- ind.1 + 1
              } else{
                ind.2 <- ind.1 - 1
              }
            }
          }
          result.mat[ii, ] <-
            (result.raw$CDE[, ind.1] + result.raw$CDE[, ind.2]) / 2
        }
      }

      ### linear&normal, parametric ###
      if (gps.method == "linear&normal" |
          gps.method == "rf&normal" |
          gps.method == "nnet&normal" | gps.method == "rf&boxcox" |
          gps.method == "nnet&boxcox" |
          gps.method == "boosting&normal" |
          gps.method == "linear&boxcox") {
        # decompose the data.eval and data.chunk
        data.eval.covariates <-
          data.eval[,!(names(data.eval) %in% c("t", "y")), drop = FALSE]
        data.eval.treatment <-
          data.eval[, (names(data.eval) %in% c("t")), drop = FALSE]
        data.chunk.covariates <-
          data.chunk[,!(names(data.chunk) %in% c("t", "y")), drop = FALSE]
        data.chunk.treatment <-
          data.chunk[, (names(data.chunk) %in% c("t")), drop = FALSE]

        ## save the gps value into the result.mat:
        for (ii in 1:dim(data.eval)[1]) {
          if (verbose == TRUE) {
            print(paste(ii / dim(data.eval)[1] * 100, "percent complete", sep = ""))
          }

          data.temp <- cbind(data.frame(t = rep(
            x = data.eval.treatment[ii, ],
            times = dim(data.chunk.covariates)[1]
          )),
          data.chunk.covariates)

          result.mat[ii,] <-
            object(data.temp)
        }
        if (verbose == TRUE) {
          cat("\n")
        }
      }

      # check whether result.mat still have NA value
      if (any(is.na(result.mat))) {
        stop("Our Evaluation Matrix as an output of function gps.series.predict is invalid!")
      }

    }

    ## return vector saving pi(T|X), used in \pi(T|X)
    result.vec <-
      rep(NA, dim(data.eval)[1]) # vector to save the final result

    ### series method ###
    if (gps.method == "series") {
      x.data.eval <-
        data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
      result.raw <-
        specSeriesCDE::predictCDE(object = object,
                   xTest = x.data.eval,
                   B = grid.length)


      for (ii in 1:dim(data.eval)[1]) {
        treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]
        ind.1 <- which.min(abs(result.raw$z - treatmentValue))
        if (ind.1 == 1) {
          ind.2 <- ind.1 + 1
        } else{
          if (ind.1 == grid.length) {
            ind.2 <- ind.1 - 1

          } else{
            if (abs(result.raw$z[ind.1 + 1] - treatmentValue) <
                abs(result.raw$z[ind.1 - 1] - treatmentValue)) {
              ind.2 <- ind.1 + 1
            } else{
              ind.2 <- ind.1 - 1
            }
          }
        }
        result.vec[ii] <-
          (result.raw$CDE[ii, ind.1] + result.raw$CDE[ii, ind.2]) / 2
      }
    }

    ### linear&normal, parametric ###
    if (gps.method == "linear&normal" |
        gps.method == "rf&normal" |
        gps.method == "nnet&normal" | gps.method == "rf&boxcox" |
        gps.method == "nnet&boxcox" |
        gps.method == "boosting&normal" |
        gps.method == "linear&boxcox") {
      result.vec <-
        object(data.eval)
    }
    # check whether result.mat still have NA value
    if (any(is.na(result.vec))) {
      stop("Our Evaluation Vector as an output of function gps.series.predict is invalid!")
    }

    # if t.grid is degenerated, then we skip this matrix-result calculation
    if (!is.null(t.grid)) {
      result.mat <-
        matrix(nrow = length(t.grid),
               ncol = dim(data.eval)[1]) # matrix to save the final result

      ### series method ###
      if (gps.method == "series") {
        x.data.eval <-
          data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
        result.raw <-
          specSeriesCDE::predictCDE(object = object,
                     xTest = x.data.eval,
                     B = grid.length)

        for (ii in 1:length(t.grid)) {
          treatmentValue <- t.grid[ii]
          ind.1 <- which.min(abs(result.raw$z - treatmentValue))
          if (ind.1 == 1) {
            ind.2 <- ind.1 + 1
          } else{
            if (ind.1 == grid.length) {
              ind.2 <- ind.1 - 1

            } else{
              if (abs(result.raw$z[ind.1 + 1] - treatmentValue) <
                  abs(result.raw$z[ind.1 - 1] - treatmentValue)) {
                ind.2 <- ind.1 + 1
              } else{
                ind.2 <- ind.1 - 1
              }
            }
          }
          result.mat[ii, ] <-
            (result.raw$CDE[, ind.1] + result.raw$CDE[, ind.2]) / 2
        }
      }

      ### linear&normal, parametric ###
      if (gps.method == "linear&normal" |
          gps.method == "rf&normal" |
          gps.method == "nnet&normal" | gps.method == "rf&boxcox" |
          gps.method == "nnet&boxcox" |
          gps.method == "boosting&normal" |
          gps.method == "linear&boxcox") {
        # decompose the data.eval and data.chunk
        data.eval.covariates <-
          data.eval[,!(names(data.eval) %in% c("t", "y")), drop = FALSE]
        data.eval.treatment <-
          data.eval[, (names(data.eval) %in% c("t")), drop = FALSE]

        for (ii in 1:length(t.grid)) {
          if (verbose == TRUE) {
            print(paste(ii / dim(data.eval)[1] * 100, "percent complete", sep = ""))
          }

          data.temp <- cbind(data.frame(t = rep(
            x = t.grid[ii],
            times = dim(data.eval)[1]
          )),
          data.eval.covariates)

          result.mat[ii,] <-
            object(data.temp)
        }
      }

      # check whether result.mat still have NA value
      if (any(is.na(result.mat))) {
        stop("Our Evaluation Matrix as an output of function gps.series.predict is invalid!")
      }
    }

    newList <- list("matrix" = result.mat, "vector" = result.vec)

    return(newList)

    ####################### abandoned codes ########################
    # for (ii in 1:dim(data.eval)[1]) {
    #   if (verbose == TRUE) {
    #     print(paste(ii / dim(data.eval)[1] * 100, "percent complete", sep = ""))
    #   }
    #   result.vec[ii] <-
    #     object(data.eval[ii,])
    # }
    # if (verbose == TRUE) {
    #   cat("\n")
    # }
  }
