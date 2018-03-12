

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
#' @param detoured if FALSE method will be done normally;
#'             otherwise a special simulation will be examined, where CDML is chosen method, CTE is chosen model,
#'             gps.method "series" is not allowed.
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
           verbose = FALSE,
           treatment.min = -4,
           treatment.max = 4,
           detoured = FALSE)
  {
    if (detoured == TRUE) {
      if (gps.method == "linear&normal" |
          gps.method == "rf&normal" |
          gps.method == "nnet&normal" | gps.method == "rf&boxcox" |
          gps.method == "nnet&boxcox" |
          gps.method == "boosting&normal" |
          gps.method == "linear&boxcox" |
          gps.method == "quantregForest" |
          gps.method == "quantregForest&normal") {
        # decompose the data.eval and data.chunk
        data.eval.covariates.treatment <-
          data.eval[,!(names(data.eval) %in% c("y")), drop = FALSE]

        res <- object(data.eval.covariates.treatment)

        return(cbind(data.eval, z = res))
      }



    }

    if (detoured == FALSE) {
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

        ### quantreg forest method ###
        if (gps.method == "quantregForest") {
          ## return matrix saving pi(T|x) for various x, used in \int_\CX \pi(T|x)dP_x
          x.data.chunk <-
            data.matrix(data.chunk[,!(names(data.chunk) %in% c("t", "y"))])
          result.raw <-
            predict(
              object,
              x.data.chunk,
              what = function(x)
                sample(x, 100, replace = TRUE)
            )

          result <-
            apply(
              X = result.raw,
              FUN = density,
              from = treatment.min,
              to = treatment.max,
              MARGIN = 1
            )
          ## save the conditional density value into the result.mat
          for (ii in 1:dim(data.eval)[1]) {
            treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]
            result.mat[ii, ] <-
              mapply(
                result,
                FUN = function(x)
                  approx(x, xout = treatmentValue)$y
              )
          }
        }

        ### quantreg forest method ###
        if (gps.method == "quantregForest&normal") {
          ## return matrix saving pi(T|x) for various x, used in \int_\CX \pi(T|x)dP_x
          x.data.chunk <-
            data.matrix(data.chunk[,!(names(data.chunk) %in% c("t", "y"))])
          result.raw <-
            predict(
              object,
              x.data.chunk,
              what = function(x)
                sample(x, 100, replace = TRUE)
            )

          res.mean <- apply(X = result.raw,
                            FUN = mean,
                            MARGIN = 1)
          res.sd <- apply(X = result.raw,
                          FUN = sd,
                          MARGIN = 1)

          ## save the conditional density value into the result.mat
          for (ii in 1:dim(data.eval)[1]) {
            treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]
            result.mat[ii, ] <-
              dnorm(
                treatmentValue,
                mean = res.mean,
                sd = res.sd,
                log = FALSE
              )
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

      ### quantreg forest method ###
      if (gps.method == "quantregForest") {
        x.data.eval <-
          data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
        result.raw <-
          predict(
            object,
            x.data.eval,
            what = function(x)
              sample(x, 100, replace = TRUE)
          )

        result <-
          apply(
            X = result.raw,
            FUN = density,
            from = treatment.min,
            to = treatment.max,
            MARGIN = 1
          )
        ## save the conditional density value into the result.mat
        for (ii in 1:dim(data.eval)[1]) {
          treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]
          result.vec[ii] <-
            approx(result[[ii]], xout = treatmentValue)$y
        }
      }


      ### quantreg forest & normal method ###
      if (gps.method == "quantregForest&normal") {
        x.data.eval <-
          data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
        result.raw <-
          predict(
            object,
            x.data.eval,
            what = function(x)
              sample(x, 100, replace = TRUE)
          )

        res.mean <- apply(X = result.raw,
                          FUN = mean,
                          MARGIN = 1)
        res.sd <- apply(X = result.raw,
                        FUN = sd,
                        MARGIN = 1)

        ## save the conditional density value into the result.mat
        for (ii in 1:dim(data.eval)[1]) {
          treatmentValue <- data.eval[ii, (names(data.eval) %in% c("t"))]
          result.vec[ii] <-
            dnorm(
              treatmentValue,
              mean = res.mean[ii],
              sd = res.sd[ii],
              log = FALSE
            )
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

        ### quantreg forest method ###
        if (gps.method == "quantregForest") {
          x.data.eval <-
            data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
          result.raw <-
            predict(
              object,
              x.data.eval,
              what = function(x)
                sample(x, 100, replace = TRUE)
            )

          result <-
            apply(
              X = result.raw,
              FUN = density,
              from = treatment.min,
              to = treatment.max,
              MARGIN = 1
            )
          ## save the conditional density value into the result.mat
          for (ii in 1:length(t.grid)) {
            treatmentValue <- t.grid[ii]
            result.mat[ii, ] <-
              mapply(
                result,
                FUN = function(x)
                  approx(x, xout = treatmentValue)$y
              )
          }
        }

        ### quantreg forest method ###
        if (gps.method == "quantregForest&normal") {
          x.data.eval <-
            data.matrix(data.eval[,!(names(data.eval) %in% c("t", "y"))])
          result.raw <-
            predict(
              object,
              x.data.eval,
              what = function(x)
                sample(x, 100, replace = TRUE)
            )

          res.mean <- apply(X = result.raw,
                            FUN = mean,
                            MARGIN = 1)
          res.sd <- apply(X = result.raw,
                          FUN = sd,
                          MARGIN = 1)

          ## save the conditional density value into the result.mat
          for (ii in 1:length(t.grid)) {
            treatmentValue <- t.grid[ii]
            result.mat[ii, ] <-
              dnorm(
                treatmentValue,
                mean = res.mean,
                sd = res.sd,
                log = FALSE
              )
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
    }
  }
