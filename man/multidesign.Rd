% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_generic.R, R/multidesign.R
\name{multidesign}
\alias{multidesign}
\alias{multidesign.matrix}
\title{Create a Multidesign Object}
\usage{
multidesign(x, y, ...)

\method{multidesign}{matrix}(x, y, column_design = NULL, ...)
}
\arguments{
\item{x}{A numeric matrix where rows are observations and columns are variables}

\item{y}{A data frame containing design variables for each observation (must have same number of rows as x)}

\item{...}{Additional arguments passed to methods, such as column_design}

\item{column_design}{Optional data frame containing metadata for columns in x (must have same number of rows as ncol(x))}
}
\value{
A multidesign object with components:
  \item{x}{The input data matrix}
  \item{design}{A tibble containing design variables}
  \item{column_design}{A tibble containing column metadata (if provided)}

A multidesign object with components:
  \item{x}{The input data matrix}
  \item{design}{A tibble containing design variables with an added .index column}
  \item{column_design}{A tibble containing column metadata (empty if not provided)}
}
\description{
Constructs a new multivariate design object linking vector-valued observations with
design variables. A multidesign object maintains the relationship between experimental
data (observations) and metadata about experimental conditions (design variables).

Creates a multidesign object that combines experimental data (as a matrix) with design information
(as a data frame) and optional column metadata. This structure is particularly useful for
experimental designs where observations have multiple associated factors and variables may have metadata.
}
\details{
A multidesign object consists of three main components:
* A data matrix where rows represent observations and columns represent variables
* A design data frame containing experimental factors and conditions for each observation
* Optional column metadata describing properties of each variable

A multidesign object consists of three main components:
* A data matrix where rows represent observations and columns represent variables
* A design data frame containing experimental factors and conditions for each observation
* Optional column metadata describing properties of each variable

The object maintains the relationship between these components while providing methods
for manipulation, subsetting, and analysis.
}
\examples{
# Create example data matrix
X <- matrix(rnorm(20*100), 20, 100)

# Create design information
Y <- tibble(condition=rep(letters[1:5], 4))

# Create multidesign object
mds <- multidesign(X, Y)

# Split by condition
sdes <- split(mds, condition)

# Create example data matrix
X <- matrix(rnorm(20*100), 20, 100)

# Create design information
Y <- tibble(
  condition = rep(c("control", "treatment"), each=10),
  subject = rep(1:5, times=4)
)

# Create column metadata
col_info <- data.frame(
  roi = paste0("region_", 1:100),
  hemisphere = rep(c("left", "right"), 50)
)

# Create multidesign object
mds <- multidesign(X, Y, col_info)

}
\seealso{
\code{\link{reduce.multidesign}} for dimensionality reduction,
  \code{\link{split.multidesign}} for splitting by design variables,
  \code{\link{multiframe}} for an alternative implementation

\code{\link{reduce.multidesign}} for dimensionality reduction,
  \code{\link{split.multidesign}} for splitting by design variables

Other multidesign functions: 
\code{\link{design.multidesign}()},
\code{\link{reduce.multidesign}()},
\code{\link{split.multidesign}()},
\code{\link{split_indices.multidesign}()},
\code{\link{subset.multidesign}()},
\code{\link{summarize_by.multidesign}()},
\code{\link{xdata.multidesign}()}
}
\concept{multidesign functions}
