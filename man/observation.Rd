% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R
\name{observation}
\alias{observation}
\title{Create an Observation Object}
\usage{
observation(x, i)
}
\arguments{
\item{x}{The data source (matrix, list, vector, or other supported object)}

\item{i}{The index of the observation to extract}
}
\value{
A function that, when called, returns the specified observation from the data source
}
\description{
Constructs a new lazy-evaluated observation object from various data sources.
An observation represents a single row, element, or vector from a data source
that can be accessed on demand.
}
\examples{
# From a matrix
X <- matrix(1:20, 5, 4)
obs1 <- observation(X, 2)  # Create observation for row 2
obs1()  # Retrieve the observation (returns row 2 as a 1xn matrix)

# From a list
X_list <- list(a=1:5, b=6:10, c=11:15)
obs2 <- observation(X_list, 3)  # Create observation for 3rd element
obs2()  # Retrieve the observation (returns the 3rd element: 11:15)

}
\seealso{
[multiframe()], [obs_group()]
}
