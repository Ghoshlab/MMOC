% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustStruct.R
\name{clustStruct}
\alias{clustStruct}
\title{Generate multi-view data sets with simple cluster structures}
\usage{
clustStruct(n, p, k, noiseDat = "random", randNoise = 2)
}
\arguments{
\item{n}{An integer, the sample size for all generated data sets}

\item{p}{An integer, the number of columns (features) in each generated data set}

\item{k}{An integer or vector, the number of distinct clusters in each generated data set. \code{n/k} must be an integer for all values of \code{k}}

\item{noiseDat}{Either the character string \code{'random'}, indicating the covariance matrix is a diagonal matrix with \code{randNoise} along the diagonal, or a valid covariance matrix}

\item{randNoise}{The value along the diagonal when \code{noiseDat='random'}}
}
\value{
A list of n\eqn{\times}p data frames with the specified number of groups
}
\description{
Generates multiple data sets from a multivariate normal distribution using the \link[MASS]{mvrnorm} function from the \code{MASS} package.
}
\details{
The function accepts \code{k} as a vector. It splits data into \code{k} groups with means \code{c(0, 2^( 1:(kk-1) ) )}, e.g., when \code{k=3} the data will be split into 3 groups with means 0, 2, and 4, respectively. The covariance matrix is either a diagonal matrix with \code{randNoise} (an integer) along the diagonal, or a given matrix.
}
\examples{

## A single view with 30 variables and 3 groups
s1 <- clustStruct(n=120, p=30, k=3, noiseDat='random')[[1]]

## Multiple views with 30 variables
## View 1 has 2 groups and View 2 has 3 groups
s2 <- clustStruct(n=120, p=30, k=c(2,3), noiseDat='random')

## Multiple views with 30 variables
## View 1 has 2 groups, View 2 has 3, and View 3 has 3 groups
s3 <- clustStruct(n=120, p=30, k=c(2,3,3), noiseDat='random')

## Three view study.
# View 1: 2 groups, 30 variables, random noise = 5
# View 2: 3 groups, 60 variables, random noise = 2
# View 3: 4 groups, 45 variables, random noise = 4

s4 <- clustStruct(n=120, k=c(2,3,4), p=c(30,60,45), randNoise=c(5,2,4))

}
