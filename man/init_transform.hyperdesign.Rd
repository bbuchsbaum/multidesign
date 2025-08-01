% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperdesign.R
\name{init_transform.hyperdesign}
\alias{init_transform.hyperdesign}
\title{Initialize Transformation for Hyperdesign}
\usage{
\method{init_transform}{hyperdesign}(x, preproc, ...)
}
\arguments{
\item{x}{A hyperdesign object}

\item{preproc}{A preprocessing specification (e.g., from recipes package)}

\item{...}{Additional arguments (not used)}
}
\value{
A list of initialized preprocessing objects, one for each block
}
\description{
Method to initialize transformations (e.g., scaling, centering) for hyperdesign objects.
Each block in the hyperdesign gets its own transformation object.
}
\seealso{
Other hyperdesign functions: 
\code{\link{design.hyperdesign}()},
\code{\link{df_to_hyperdesign}()},
\code{\link{hyperdesign}()},
\code{\link{subset.hyperdesign}()},
\code{\link{xdata.hyperdesign}()}
}
\concept{hyperdesign functions}
