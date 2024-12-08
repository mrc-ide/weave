% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel.R
\name{periodic_kernel}
\alias{periodic_kernel}
\title{Modified periodic kernel with long-term decay}
\usage{
periodic_kernel(distance, periodic_scale, long_term_scale, period)
}
\arguments{
\item{periodic_scale}{The length scale parameter for the periodic kernel (in time units)}

\item{long_term_scale}{The length scale parameter for the long term kernel (in time units)}

\item{period}{The period of the data, representing the time interval over which periodic behavior repeats (in time units)}

\item{t}{A scalar representing the time point(s) for the first input}

\item{t_prime}{A scalar or vector representing the time point(s) for the second input}
}
\description{
Period kernel refs https://tinyurl.com/52mzsrhh, https://tinyurl.com/yc4dtm68
}
