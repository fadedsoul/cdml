#' Plot of Alternatives
#'
#' @param gps.method a vector of method for generalized propensity score estimation
#' @param g.method a vector of method for regression estimation
#' @param method.involved e.g.: c("SR", "CDML.1", "CDML.2", "CDML.5", "HI")
#' @param sample.size.option a vector of number of samples that will be selected from
#' @param num.cov.option a vector of number of covariates that will be selected from
#' @param treatment.sd.option a vector of standard deviation of treatment generation that will be selected from
#' @param struct.model.option a vector of parametric curves of response generation, typically choose from "linear", "polynom", "polynom2", "polynom3", "mixture".
#' @param what choose from "sample.sizes", "cov"
#' @param file file that we save our return into
#' @param main if main == TRUE, we plot a new one, otherwise we attach our figure in an existing plot
#' @param plot.what choose from "bias" and "error"
#' @param color can be used to specify the color of the main line
#'
#' @return plots
#' @export
#'
plot.alternative.inside <-
  function(gps.method,
           g.method,
           # method = c("SR", "CDML.1", "CDML.2", "CDML.5", "HI"),
           method.involved = c("CDML.1", "CDML.2"),
           sample.size.option = NULL,
           num.cov.option = c(2, 5),
           treatment.sd.option = NULL,
           struct.model.option = NULL,
           what,
           file = NULL,
           plot.what = "bias",
           main = TRUE,
           color = "red") {
    ############# pre-processing #############

    res <- readRDS(file = file)

    error.list.list <- res$error.list.list
    bias.list <- res$bias.list


    for (i in 1:length(error.list.list)) {
      error.list.list[[i]]$error1 <-
        round(mean(
          as.numeric(error.list.list[[i]]$error1),
          na.rm = TRUE,
          trim = 0
        ), 2)
      error.list.list[[i]]$error2 <-
        round(mean(
          as.numeric(error.list.list[[i]]$error2),
          na.rm = TRUE,
          trim = 0
        ), 2)
    }

    ########### different cases ################
    if (what == "cov") {
      result.tableout.error1 <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(num.cov.option))
      result.tableout.error2 <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(num.cov.option))
      result.tableout.bias <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(num.cov.option))

      for (tt in 1:length(method.involved)) {
        ### select the indices that correspondent to our interest ###
        for (ss in 1:length(num.cov.option)) {
          array <- NULL

          if (method.involved[tt] == "CDML.1") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.2") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 2)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.3") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 3)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.5") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 5)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "SR") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$method == "SR" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "HI") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$covariates == num.cov.option[ss] &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "HI" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }
          ### the average value of error1,2 respectively,
          ## and then put it into the result.table
          for (ii in 1:length(array)) {
            result.tableout.error1[tt, ss]  <-
              1 / length(array) * error.list.list[[array[ii]]]$error1 + result.tableout.error1[tt, ss]
            result.tableout.error2[tt, ss]  <-
              1 / length(array) * error.list.list[[array[ii]]]$error2 + result.tableout.error2[tt, ss]
            result.tableout.bias[tt, ss]  <-
              1 / length(array) *  bias.list[[array[ii]]]$bias + result.tableout.bias[tt, ss]
          }
        }
      }


      ################ plotting part #################

      #### plot of errors ####
      if(plot.what == "error"){

      ylimit.error <- c(min(
        min(result.tableout.error1),
        min(result.tableout.error2)
      ), min(max(
        max(result.tableout.error1),
        max(result.tableout.error2)
      ), 2000))

      color.array <-
        c(color, "green", "blue", "purple", "cyan", "black")

      if(main == TRUE){
      plot(
        result.tableout.error1[1, ] ~ num.cov.option,
        type = "b" ,
        bty = "l" ,
        xlab = "number of covariates" ,
        ylab = "error" ,
        col = "red" ,
        lwd = 2 ,
        pch = 15,
        ylim = ylimit.error,
        log = "y"
        #xaxt = "n"
      )
}
      for (rr in 1:length(method.involved)) {
        lines(
          result.tableout.error1[rr,] ~ num.cov.option,
          col = color.array[rr] ,
          lwd = 3 ,
          pch = 19 ,
          type = "b"
        )

        points(
          result.tableout.error2[rr, ] ~ num.cov.option,
          col = color.array[rr],
          lwd = 2 ,
          pch = 16 ,
          lty = "dotted",
          type = "b"
        )
      }
}

      #### plot of bias ####
      if(plot.what == "bias"){
      ylimit.bias <- c(min(result.tableout.bias)
                       , max(result.tableout.bias))

      color.array <-
        c(color, "green", "blue", "purple", "cyan", "black")

      if(main == TRUE){
      plot(
        result.tableout.bias[1, ] ~ num.cov.option ,
        type = "b" ,
        bty = "l" ,
        xlab = "number of covariates" ,
        ylab = "bias" ,
        col = "red" ,
        lwd = 2 ,
        pch = 15,
        ylim = ylimit.bias,
        log = "y"
        #xaxt = "n"
      )
      }
      lines(
        result.tableout.bias[1,] ~ num.cov.option,
        col = color.array[1] ,
        lwd = 3 ,
        pch = 19 ,
        type = "b"
      )
}
    }

    if (what == "sample.sizes") {
      result.tableout.error1 <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(sample.size.option))
      result.tableout.error2 <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(sample.size.option))
      result.tableout.bias <-
        matrix(0,
               nrow = length(method.involved),
               ncol = length(sample.size.option))

      for (tt in 1:length(method.involved)) {
        ### select the indices that correspondent to our interest ###
        for (ss in 1:length(sample.size.option)) {
          array <- NULL

          if (method.involved[tt] == "CDML.1") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.2") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 2)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.3") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 3)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "CDML.5") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "CDML" &&
                  error.list.list[[i]]$record$fold == 5)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "SR") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$g.method == g.method &&
                  error.list.list[[i]]$record$method == "SR" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }

          if (method.involved[tt] == "HI") {
            for (i in 1:length(error.list.list)) {
              if (error.list.list[[i]]$record$samples == sample.size.option[ss] &&
                  error.list.list[[i]]$record$gps.method == gps.method &&
                  error.list.list[[i]]$record$method == "HI" &&
                  error.list.list[[i]]$record$fold == 1)
                array <- c(array, i)
            }
          }

          ### the average value of error1,2 respectively,
          ## and then put it into the result.table
          for (ii in 1:length(array)) {
            result.tableout.error1[tt, ss]  <-
              1 / length(array) * error.list.list[[array[ii]]]$error1 + result.tableout.error1[tt, ss]
            result.tableout.error2[tt, ss]  <-
              1 / length(array) * error.list.list[[array[ii]]]$error2 + result.tableout.error2[tt, ss]
            result.tableout.bias[tt, ss]  <-
              1 / length(array) *  bias.list[[array[ii]]]$bias + result.tableout.bias[tt, ss]
          }
        }
      }


      ################ plotting part #################

      if (plot.what == "error") {
        #### plot of errors ####
        ylimit.error <- c(min(
          min(result.tableout.error1),
          min(result.tableout.error2)
        ), min(max(
          max(result.tableout.error1),
          max(result.tableout.error2)
        ), 2000))

        color.array <-
          c(color, "green", "blue", "purple", "cyan", "black")

        if (main == TRUE) {
          plot(
            result.tableout.error1[1, ] ~ sample.size.option,
            type = "b" ,
            bty = "l" ,
            xlab = "number of samples" ,
            ylab = "error" ,
            col = "red" ,
            lwd = 2 ,
            pch = 15,
            ylim = ylimit.error,
            log = "y"
            #xaxt = "n"
          )
        }
        for (rr in 1:length(method.involved)) {
          lines(
            result.tableout.error1[rr,] ~ sample.size.option,
            col = color.array[rr] ,
            lwd = 3 ,
            pch = 19 ,
            type = "b"
          )

          points(
            result.tableout.error2[rr, ] ~ sample.size.option,
            col = color.array[rr],
            lwd = 2 ,
            pch = 16 ,
            lty = "dotted",
            type = "b"
          )
        }
      }

      #### plot of bias ####
      if (plot.what == "bias") {
        ylimit.bias <- c(min(result.tableout.bias)
                         , max(result.tableout.bias))

        color.array <-
          c(color, "green", "blue", "purple", "cyan", "black")


        if (main == TRUE) {
          plot(
            result.tableout.bias[1, ] ~ sample.size.option ,
            type = "b" ,
            bty = "l" ,
            xlab = "number of samples" ,
            ylab = "bias" ,
            col = "red" ,
            lwd = 2 ,
            pch = 15,
            ylim = ylimit.bias,
            log = "y"
            #xaxt = "n"
          )
        }

        lines(
          result.tableout.bias[1,] ~ sample.size.option,
          col = color.array[1] ,
          lwd = 3 ,
          pch = 19 ,
          type = "b"
        )
      }
    }

  }
