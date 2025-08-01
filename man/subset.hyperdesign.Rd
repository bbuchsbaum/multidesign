% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperdesign.R
\name{subset.hyperdesign}
\alias{subset.hyperdesign}
\title{Subset a Hyperdesign Object}
\usage{
\method{subset}{hyperdesign}(x, fexpr)
}
\arguments{
\item{x}{A hyperdesign object}

\item{fexpr}{Filter expression to apply to each block's design}
}
\value{
A new hyperdesign object containing only the selected blocks
}
\description{
Create a new hyperdesign object containing only the specified blocks.
}
\examples{
# Create example hyperdesign
d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=1, condition=rep(c("A","B"), 5)))
d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=2, condition=rep(c("A","B"), 5)))
d3 <- multidesign(matrix(rnorm(10*20), 10, 20),
                  data.frame(subject=3, condition=rep(c("A","B"), 5)))
hd <- hyperdesign(list(d1, d2, d3))

# Keep only blocks 1 and 3
subset_hd <- subset(hd, c(TRUE, FALSE, TRUE))
}
\seealso{
Other hyperdesign functions: 
\code{\link{design.hyperdesign}()},
\code{\link{df_to_hyperdesign}()},
\code{\link{hyperdesign}()},
\code{\link{init_transform.hyperdesign}()},
\code{\link{xdata.hyperdesign}()}
}
\concept{hyperdesign functions}
