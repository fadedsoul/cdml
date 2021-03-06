% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gps.estimation.R
\name{gps.method.estimation}
\alias{gps.method.estimation}
\title{Generalized Propensity Score Estimation}
\usage{
gps.method.estimation(data, gps.method = "series", n.prop = 0.7,
  bw.length.out = 10, bw.from = 0.5, bw.to = 9.5, treatment.min = -10,
  treatment.max = 10, verbose = TRUE, detoured = FALSE)
}
\arguments{
\item{data}{data}

\item{gps.method}{methods for generalized propensity score estimation, choose from "series", "linear&normal", "boxcox"}

\item{n.prop}{default = 0.7 meaning 70\% for training usage, 30\% for validation usage in method "series"}

\item{bw.length.out}{number of tested points of tuning parameters from the kernel (bandwidth) in method "series"}

\item{bw.from}{lower bound of tuning parameters from the kernel (bandwidth) in method "series"}

\item{bw.to}{upper bound of tuning parameters from the kernel (bandwidth) in method "series"}

\item{treatment.min}{in what interval conditional density shall be evaluated}

\item{treatment.max}{in what interval conditional density shall be evaluated}

\item{verbose}{print in console the procedure step by step}

\item{detoured}{if FALSE method will be done normally;
otherwise a special simulation will be examined, where CDML is chosen method, CTE is chosen model,
gps.method "series" is not allowed.}
}
\value{
object of generalized propensity score estimation
}
\description{
Generalized Propensity Score Estimation
}
