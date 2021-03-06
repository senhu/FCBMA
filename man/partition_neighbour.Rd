% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition_neighbour.R
\name{partition_neighbour}
\alias{partition_neighbour}
\alias{partition_random_neighbour}
\alias{partition_all_neighbour}
\alias{partition_all_neighbour_restricted}
\title{Generate neighbours of a graycode with or without pre-specified restrictions}
\usage{
partition_random_neighbour(x, method = "ChangeOne",
  IncludeOrigin = FALSE)

partition_all_neighbour(x, method = "ChangeOne", IncludeOrigin = FALSE)

partition_all_neighbour_restricted(x, group = NULL, apart = NULL,
  method = "ChangeOne", IncludeOrigin = FALSE)
}
\arguments{
\item{x}{a numerical vector representing a graycode}

\item{method}{which method used to generate neighbours; either \code{"ChangeOne"} or \code{"GroupSplit"}}

\item{IncludeOrigin}{logical; indicating whether the input graycode should be included}

\item{group}{a vector indicating which elements should be grouped together}

\item{apart}{a vector indicating which elements should not be grouped together}
}
\value{
a random neighbouring partition or a matrix of all neighbouring partitions
}
\description{
A graycode neighbour is generated by randomly picking two groups and put them together, or randomly pick one group and split it into 2 groups, only generate one random neighbour each time
}
\examples{
partition_random_neighbour(x = c(1,2,3,1),
                           IncludeOrigin = FALSE) # "ChangeOne" method by default
partition_all_neighbour(x = c(1,2,3,1),
                        IncludeOrigin = FALSE)
partition_random_neighbour(c(1,2,3,4,2,4,1),
                           method="GroupSplit",
                           IncludeOrigin = FALSE)
partition_all_neighbour(c(1,2,2,3,3,3,3),
                        method="GroupSplit",
                        IncludeOrigin = FALSE)
partition_all_neighbour_restricted(x = c(1,2,3,1),
                                   group=NULL,
                                   apart=c(2,3),
                                   IncludeOrigin = TRUE)
partition_all_neighbour_restricted(x=c(1,2,3,1),
                                   group = c(2,3),
                                   apart = NULL,
                                   IncludeOrigin = TRUE)
\donttest{
# the following example shows the distribution of generating random partition is not uniform
x <- c(1,2,2,2,3,3,3)
allneighbour <- NULL
for (j in c(1:1000)){
  oneneighbour <- partition_random_neighbour(x, method = "GroupSplit")
  allneighbour <- rbind(allneighbour, oneneighbour)
}
dim(allneighbour)
uniquecode <- unique(allneighbour)
dim(uniquecode)
dim(partition_all_neighbour(x, method="GroupSplit"))
count1<-count2<-count3<-count4<-count5<-count6<-count7<-count8<-count9<-0
for (i in c(1: dim(allneighbour)[1])){
  if (all(uniquecode[1,] == allneighbour[i,])) count1 <- count1 + 1
  if (all(uniquecode[2,] == allneighbour[i,])) count2 <- count2 + 1
  if (all(uniquecode[3,] == allneighbour[i,])) count3 <- count3 + 1
  if (all(uniquecode[4,] == allneighbour[i,])) count4 <- count4 + 1
  if (all(uniquecode[5,] == allneighbour[i,])) count5 <- count5 + 1
  if (all(uniquecode[6,] == allneighbour[i,])) count6 <- count6 + 1
  if (all(uniquecode[7,] == allneighbour[i,])) count7 <- count7 + 1
  if (all(uniquecode[8,] == allneighbour[i,])) count8 <- count8 + 1
  if (all(uniquecode[9,] == allneighbour[i,])) count9 <- count9 + 1
}
c(count1, count2, count3, count4, count5, count6, count7, count8, count9)
}

}
