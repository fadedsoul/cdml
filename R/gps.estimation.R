


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
#' @param detoured if FALSE method will be done normally;
#'             otherwise a special simulation will be examined, where CDML is chosen method, CTE is chosen model,
#'             gps.method "series" is not allowed.
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
                                  verbose = TRUE,
                                  detoured = FALSE)

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
      data[indices, ][, !(names(data) %in% c("y"))] # training data
    data.valid <-
      data[-indices, ][, !(names(data) %in% c("y"))] # validation data

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
          xTrain = data.matrix(data.train[, !(names(data.train) %in% c("t"))]),
          xValidation = data.matrix(data.valid[, !(names(data.valid) %in% c("t"))]),
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
        xTrain = data.matrix(data.train[, !(names(data.train) %in% c("t"))]),
        xValidation = data.matrix(data.valid[, !(names(data.valid) %in% c("t"))]),
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
    data2 <- data[, !(names(data) %in% c("y"))]
    lm.res <- lm(t ~ ., data = data2)

    mean <- mean(lm.res$residuals)
    sd <- sd(lm.res$residuals)

    if(detoured == FALSE){
    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(lm.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
    }

    if(detoured == TRUE){
      # data <- cbind( data, z = lm.res$residuals)
      object <- function(data.eval) {
        return(
          data.eval$t - predict(lm.res, data.eval)
          )
      }
    }
  }

  #### rf model of t,x, and the noise is normal distributed ####
  if (gps.method == "rf&normal") {
    # filter out the y column
    data2 <- data[, !(names(data) %in% c("y"))]
    rf.res <- randomForest(t ~ ., data = data2, ntrees = 1000)

    mean <- mean(data$t - rf.res$predicted)
    sd <- sd(data$t - rf.res$predicted)

    if(detoured == FALSE){
    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(rf.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
    }

    if(detoured == TRUE){
      data <- cbind( data, z = (data$t - rf.res$predicted))
    }
  }

  #### nnet model of t,x, and the noise is normal distributed ####
  if (gps.method == "nnet&normal") {
    # filter out the y column
    data2 <- data[, !(names(data) %in% c("y"))]
    nnet.res <- nnet(t ~ .,
                     size = 8,
                     linout = TRUE,
                     data = data2)

    mean <- mean(nnet.res$residuals)
    sd <- sd(nnet.res$residuals)

    if(detoured == FALSE){
    object <- function(data.eval) {
      return(dnorm(
        data.eval$t - predict(nnet.res, data.eval),
        mean = mean,
        sd = sd,
        log = FALSE
      ))
    }
    }

    if(detoured == TRUE){
      data <- cbind( data, z = (data$t - nnet.res$predicted))
    }

  }

  #### boosting model of t,x and the noise is normal distributed ####
  if (gps.method == "boosting&normal") {
    # filter out the y column
    data <- data[, !(names(data) %in% c("y"))]

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
        newsample = data[bo,]
        j.drop = match(c("t"), names(data))
        j.drop = j.drop[!is.na(j.drop)]
        x = newsample[,-j.drop]
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

    res <-
      data$t - predict(
        model.den,
        newdata = data,
        n.trees = floor(best.aac.iter),
        type = "response"
      )
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
                   mean = mean, sd = sd))
    }
  }

  #### quantile regression Random Forest method ####
  ### www.jmlr.org/papers/v7/meinshausen06a.html ####
  if(gps.method == "quantregForest&normal"){
    Xtrain     <- data[, !(names(data) %in% c("t", "y"))]
    Ytrain     <- data[, names(data) %in% c("t")]

    if(detoured == FALSE){
      object <- quantregForest(x = Xtrain, y = Ytrain)
    }

    if(detoured == TRUE){

      object2 <- quantregForest(x = Xtrain, y = Ytrain)


      object <- function(data.eval) {

        x.data <-
          data.matrix(data.eval[, !(names(data.eval) %in% c("y"))])

        result.raw <-
          predict(object2, x.data,what = function(x) sample(x,100,replace=TRUE))

        result <- apply(X=result.raw, FUN = mean, MARGIN = 1)

        return(
          data.eval$t - as.vector(result)
        )

      }
    }


    }

  ####neural network box-cox ####
  if (gps.method == "nnet&boxcox") {
    ## filter out the y column
    data <- data[, !(names(data) %in% c("y"))]
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
    data <- data[, !(names(data) %in% c("y"))]
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
    data <- data[, !(names(data) %in% c("y"))]
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

  #### quantile regression Random Forest method ####
  ### www.jmlr.org/papers/v7/meinshausen06a.html ####
  if (gps.method == "quantregForest") {
    Xtrain     <- data[, !(names(data) %in% c("t", "y"))]
    Ytrain     <- data[, names(data) %in% c("t")]

    if(detoured == FALSE){
    object <- quantregForest(x = Xtrain, y = Ytrain)
    }
    if(detoured == TRUE){

      object2 <- quantregForest(x = Xtrain, y = Ytrain)


      object <- function(data.eval) {

        x.data <-
          data.matrix(data.eval[, !(names(data.eval) %in% c("y"))])

        result.raw <-
          predict(object2, x.data,what = function(x) sample(x,100,replace=TRUE))

        result <- apply(X=result.raw, FUN = mean, MARGIN = 1)

        return(
          data.eval$t - as.vector(result)
        )

      }
    }
  }

  ## what to return finally
  return(object)
}
