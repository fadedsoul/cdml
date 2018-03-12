plot.alternative.inside <-
  function(gps.method,
           g.method,
           # method = c("SR", "CDML.1", "CDML.2", "CDML.5", "HI"),
           method.involved = c("CDML.1","CDML.2"),
           sample.size.option = NULL,
           num.cov.option = c(2, 5),
           treatment.sd.option = NULL,
           struct.model.option = NULL,
           what,
           file = NULL) {
    ############# pre-processing #############

    try(load(file = file))

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
        matrix(0, nrow = length(method.involved), ncol = length(num.cov.option))
      result.tableout.error2 <-
        matrix(0, nrow = length(method.involved), ncol = length(num.cov.option))
      result.tableout.bias <-
        matrix(0, nrow = length(method.involved), ncol = length(num.cov.option))

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

    }
    ################ plotting part #################

    #### plot of errors ####
    ylimit.error <- c(min(min(result.tableout.error1),
                          min(result.tableout.error2)), min(max(
                            max(result.tableout.error1),
                            max(result.tableout.error2)
                          ), 2000))

    color.array <- c("red", "green", "blue", "purple", "cyan", "black")

    plot(
      result.tableout.error1[1,] ~ num.cov.option,
      type = "b" ,
      bty = "l" ,
      xlab = "covariance" ,
      ylab = "" ,
      col = "red" ,
      lwd = 2 ,
      pch = 15,
      ylim = ylimit.error,
      log = "y",
      xaxt = "n"
    )

    for(rr in 1:length(method.involved)){

      lines(
        result.tableout.error1[rr, ] ~ num.cov.option,
        col = color.array[rr] ,
        lwd = 3 ,
        pch = 19 ,
        type = "b"
      )

      points(
        result.tableout.error2[rr,] ~ num.cov.option,
        col = color.array[rr],
        lwd = 2 ,
        pch = 16 ,
        lty = "dotted",
        type = "b"
      )
    }


    #### plot of bias ####
    ylimit.bias <- c(min(result.tableout.bias)
                     , max(result.tableout.bias))

    plot(
      result.tableout.bias[1,] ~ num.cov.option ,
      type = "b" ,
      bty = "l" ,
      xlab = "se" ,
      ylab = "" ,
      col = "red" ,
      lwd = 2 ,
      pch = 15,
      ylim = ylimit.bias,
      log = "y",
      xaxt = "n"
    )
  }
