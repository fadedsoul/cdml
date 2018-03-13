setwd("/sfs/u/zgong/RsourceCode")

plot.what = "bias"

num.cov.option = c(2, 3, 5, 10)
sample.size.option = c(200, 300,500, 750, 1000, 1500)

plot.alternative.inside(
  gps.method = "quantregForest&normal",
  g.method = "rf",
  num.cov.option = c(2, 3, 5, 10),
  what = "cov",
  method.involved = c("CDML.1"),
  file = "qrf_1_detoured.rds",
  plot.what = plot.what
)
plot.alternative.inside(
  gps.method = "quantregForest&normal",
  g.method = "rf",
  num.cov.option = c(2, 3, 5, 10),
  what = "cov",
  method.involved = c("CDML.1"),
  file = "qrf_1.rds",
  plot.what = plot.what,
  main = FALSE,
  color = "black"
)

plot.alternative.inside(
  gps.method = "rf&normal",
  g.method = "rf",
  num.cov.option = c(2, 3, 5, 10),
  what = "cov",
  method.involved = c("CDML.1"),
  file = "rf_1_detoured.rds",
  plot.what = plot.what,
  main = FALSE,
  color = "blue"
)

plot.alternative.inside(
  gps.method = "rf&normal",
  g.method = "rf",
  num.cov.option = c(2, 3, 5, 10),
  what = "cov",
  method.involved = c("CDML.1"),
  file = "rf_1.rds",
  plot.what = plot.what,
  main = FALSE,
  color = "cyan"
)

plot.alternative.inside(
  gps.method = "quantregForest&normal",
  g.method = "rf",
  num.cov.option = c(2, 3, 5, 10),
  what = "cov",
  method.involved = c("CDML.2"),
  file = "qrf_2_detoured.rds",
  plot.what = plot.what,
  main = FALSE,
  color = "green"
)
