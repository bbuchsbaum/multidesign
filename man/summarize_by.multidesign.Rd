% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multidesign.R
\name{summarize_by.multidesign}
\alias{summarize_by.multidesign}
\title{Summarize a Multidesign Object by Design Variables}
\usage{
\method{summarize_by}{multidesign}(x, ..., sfun = colMeans, extract_data = FALSE)
}
\arguments{
\item{x}{A multidesign object}

\item{...}{Unquoted names of variables to group by}

\item{sfun}{Summary function to apply (default is colMeans)}

\item{extract_data}{Logical; whether to extract raw data instead of computing summary}
}
\value{
A new multidesign object containing:
  \item{x}{Matrix of summary statistics}
  \item{design}{Design information for each summary}
  \item{column_design}{Original column metadata}
}
\description{
Computes summaries of the data matrix grouped by combinations of design variables.
}
\examples{
X <- matrix(rnorm(100*20), 100, 20)
Y <- tibble(
  condition = rep(c("A", "B"), each=50),
  block = rep(1:2, times=50)
)
mds <- multidesign(X, Y)

# Get means by condition
means_by_cond <- summarize_by(mds, condition)

# Get means by condition and block
means_by_both <- summarize_by(mds, condition, block)

}
\seealso{
\code{\link{split.multidesign}}

Other multidesign functions: 
\code{\link{design.multidesign}()},
\code{\link{multidesign}()},
\code{\link{reduce.multidesign}()},
\code{\link{split.multidesign}()},
\code{\link{split_indices.multidesign}()},
\code{\link{subset.multidesign}()},
\code{\link{xdata.multidesign}()}
}
\concept{multidesign functions}
