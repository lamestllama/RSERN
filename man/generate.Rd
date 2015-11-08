\name{generate}
\alias{generate}
\title{Spatially Embedded Random Network Generator}
\description{Generates a Spatially Embedded Random Network}
\usage{
generate(s, q, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{s parameter}
  \item{q}{q parameter}
  \item{N}{number of nodes}
}
\value{
An igraph object including
\item{weight}{length of each edge}
\item{x}{x coordinates of nodes}
\item{y}{y coordinates of nodes}
}
\references{
% TODO Add reference
}
\author{Eric Parsonage}
\examples{
x <- generate(s=1,q=0.1,N=1000)
}
