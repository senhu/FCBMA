% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graycode_to_partition.R
\name{graycode_to_partition}
\alias{graycode_to_partition}
\title{Convert graycode to listing of partition}
\usage{
graycode_to_partition(x)
}
\arguments{
\item{x}{a graycode (in any format, not necessarily in canonical format)}
}
\value{
a character string listing equivalence class
}
\description{
Convert a graycode to listing of partition (equivalence class), to illustrate which elements are grouped together.
}
\examples{
graycode_to_partition(c(1,2,3,1,2))
graycode_to_partition(c(2,1,3,2,1))

}
\seealso{
\code{\link[partitions]{listParts}}
}
