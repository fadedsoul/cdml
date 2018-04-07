

################################################################################
#                                                                              #
#                          g.method.estimation                                 #
#                                                                              #
################################################################################

#' Regression Method Estimation
#'
#' @param data data that we based estimation on.
#' @param g.method choose regression method from "nnet", "trees", "rf", "MARS", "lasso".
#'
#' @return results of correspondent method estimation
#'
g.method.estimation <- function(data, g.method) {
  ### random forest ###
  if (g.method == "rf") {
    g.rf <- randomForest(y ~ ., data = data, ntrees = 500)
    return(g.rf)
  }

  ### neural net ###
  if (g.method == "nnet") {
    g.nnet <- nnet(y ~ .,
                   data = data,
                   size = 2,
                   decay=0.1,
                   linout = TRUE)
    return(g.nnet)
  }

  ### trees (rpart) ###
  if (g.method == "trees") {
    trees          <-
      rpart(
        y ~ .,
        data = data,
        method = "anova",
        control = rpart.control(minsplit = 10, cp = 0.0001)
      )
    bestcp         <-
      trees$cptable[which.min(trees$cptable[, "xerror"]), "CP"]
    g.trees          <- prune(trees, cp = bestcp)
    return(g.trees)
  }

  ## MARS ###
  if (g.method == "MARS") {
    g.mars <- earth(y ~ ., data = data, thresh = 0.01, degree = 1)
    return(g.mars)
  }

  ## lasso ###
  if (g.method == "lasso") {
    g.lasso <-
      glmnet(
        y = as.matrix(data[, (names(data) %in% c("y"))]),
        x = as.matrix(data[, !(names(data) %in% c("y"))]),
        alpha = 1,
        intercept = TRUE
      )
    return(g.lasso)
  }
}

################################################################################
#                                                                              #
#                                g.predict                                     #
#                                                                              #
################################################################################

#' Regression Method Prediction
#'
#' @param object object that is from g.method.estimation function
#' @param data.eval data to be evaluated, the (T,X) pair is of interest
#' @param data.chunk whole data, covariates are of interest
#' @param mode choose from "CDML", "SR", "HI". The default set-up is for "CDML" case.
#' @param verbose whether to print step in console
#' @param g.method regression g.method that is specified in parameter object.
#'
#' @return a vector and a matrix
#' For the matrix elements, they stand for the predicted values of:
#' row means given certain treatment from data.eval
#' column means given certain covariates from data.chunk.
#' For the vector: we list in order the predicted value as \eqn{i=1,...,n} for \eqn{(T_i, X_i)}
g.predict <-
  function(object,
           data.eval,
           data.chunk,
           mode = "CDML",
           verbose = TRUE,
           g.method = "rf")
  {
    # the default set-up, we invoke this function in CDML estimation
    if (mode == "CDML") {
      # decompose the data.eval and data.chunk
      data.eval.covariates <-
        data.eval[, !(names(data.eval) %in% c("t", "y")), drop = FALSE]
      data.eval.treatment <-
        data.eval[, (names(data.eval) %in% c("t")), drop = FALSE]
      data.chunk.covariates <-
        data.chunk[, !(names(data.chunk) %in% c("t", "y")), drop = FALSE]
      data.chunk.treatment <-
        data.chunk[, (names(data.chunk) %in% c("t")), drop = FALSE]

      ## debug mode

     # print(data.eval.covariates)
     # print(data.eval.treatment)


      print(data.eval[, !(names(data.eval) %in% c("y")), drop = FALSE])
      ##



      ## return matrix saving g(T|x) for various x, used in \int_\CX g(T|x)dP_x
      result.mat <-
        matrix(nrow = dim(data.eval)[1],
               ncol = dim(data.chunk)[1]) # matrix to save the final result
      # row means given certain treatment from data.eval
      # column means given certain covariates from data.complete

      ## save the g value into the result.mat:
      for (ii in 1:dim(data.eval)[1]) {
        if (verbose == TRUE) {
          print(paste(ii / dim(data.eval)[1] * 100, "percent complete", sep = ""))
        }

        ## lasso method
        if (g.method == "lasso") {
          mat <- as.matrix(cbind(data.frame(
            t = rep(
              x = data.eval.treatment[ii,],
              times = dim(data.chunk.covariates)[1]
            )
          ),
          data.chunk.covariates))

          result.mat[ii, ] <-
            predict(object, mat, s = 0.01)

        } else{
          result.mat[ii, ] <-
            predict(object, cbind(data.frame(
              t = rep(
                x = data.eval.treatment[ii,],
                times = dim(data.chunk.covariates)[1]
              )
            ),
            data.chunk.covariates))
        }
        if (verbose == TRUE) {
          cat("\n")
        }
      }

      # check whether result.mat still have NA value
      if (any(is.na(result.mat))) {
        stop("Our Evaluation Matrix as an output of function gps.series.predict is invalid!")
      }

      ## return vector saving g(T|X), used in g(T|X)
      result.vec <-
        rep(NA, dim(data.eval)[1]) # vector to save the final result

      for (ii in 1:dim(data.eval)[1]) {
        if (verbose == TRUE) {
          print(paste(ii / dim(data.eval)[1] * 100, "percent complete", sep = ""))
        }

        ## g.method lasso
        if (g.method == "lasso") {
          mat <- as.matrix(cbind(
            data.frame(t = data.eval.treatment[ii,]),
            data.eval.covariates[ii,]
          ))

          result.vec[ii] <-
            predict(object, mat, s = 0.01)

        } else{
          # result.vec[ii] <-
          #   predict(object, cbind(
          #     data.frame(t = data.eval.treatment[ii,]),
          #     data.eval.covariates[ii,]
          #   ))

          result.vec[ii] <-
            predict(object,  data.eval[ii, !(names(data.eval) %in% c("y")), drop = FALSE])


        }
        if (verbose == TRUE) {
          cat("\n")
        }
      }

      # check whether result.mat still have NA value
      if (any(is.na(result.vec))) {
        stop("Our Evaluation Vector as an output of function gps.series.predict is invalid!")
      }

      newList <- list("matrix" = result.mat, "vector" = result.vec)

      return(newList)
    }

    if (mode == "SR") {
      # decompose the data.chunk
      data.chunk.covariates <-
        data.chunk[, !(names(data.chunk) %in% c("t", "y")), drop = FALSE]
      data.chunk.treatment <-
        data.chunk[, (names(data.chunk) %in% c("t")), drop = FALSE]

      ## return matrix saving g(T|x) for various x, used in \int_\CX g(T|x)dP_x
      result.mat <-
        matrix(nrow = length(data.eval),
               ncol = dim(data.chunk)[1]) # matrix to save the final result
      # row means given certain treatment from data.eval
      # column means given certain covariates from data.complete

      ## save the conditional density value into the result.mat:
      for (ii in 1:length(data.eval)) {
        #print(paste(ii / length(data.eval) * 100, "percent complete", sep = ""))

        if (g.method == "lasso") {
          mat <- as.matrix(cbind(data.frame(t = rep(
            x = data.eval[ii],
            times = dim(data.chunk.covariates)[1]
          )),
          data.chunk.covariates))

          result.mat[ii, ] <-
            predict(object, mat, s = 0.01)
        } else{
          result.mat[ii, ] <-
            predict(object, cbind(data.frame(t = rep(
              x = data.eval[ii],
              times = dim(data.chunk.covariates)[1]
            )),
            data.chunk.covariates))
        }
      }

      # check whether result.mat still have NA value
      if (any(is.na(result.mat))) {
        stop("Our Evaluation Matrix as an output of function g.rf.predict is invalid!")
      }

      newList <- list("matrix" = result.mat)

      return(newList)
    }
  }
