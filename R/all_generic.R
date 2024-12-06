#' create an observation object
#'
#' construct a new multivariate observation vector
#'
#' @param x the data source
#' @param i the index of the observation
#' @export
observation <- function(x, i) UseMethod("observation")



#' create an multidesign object
#'
#' construct a new multivariate design object linking vector-valued 
#' observations and design variables
#'
#' @param x the multivariate data (a matrix, a list, or other data container)
#' @param y a design matrix with same number of rows/elements as x
#' @param ... extra args
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble(condition=rep(letters[1:5], 4))
#'
#' mds <- multidesign(X,Y)
#' sdes <- split(mds, condition)
#' @export
multidesign <- function(x, y, ...) UseMethod("multidesign")

#' Create a Hyperdesign Object
#'
#' Construct a new `hyperdesign` object that encapsulates a collection of `multidesign` objects.
#' Used to model multi-group, multi-block, or multi-view datasets, where each block/group/view 
#' is associated with a matrix-variate response and an arbitrary design.
#'
#' @param x the input data (typically a list of multidesign objects)
#' @param ... additional arguments passed to methods
#' @return A `hyperdesign` object
#' @seealso [multidesign()], [multiblock()]
#' @export
hyperdesign <- function(x, ...) UseMethod("hyperdesign")



#' Create a Multiframe Object
#'
#' Construct a new multivariate design object linking vector-valued observations with 
#' arbitrary design variables. Similar to `multidesign` but with enhanced support for 
#' data frame operations.
#'
#' @param x the multivariate data (a matrix, list, or other data container)
#' @param y a design matrix or data frame with same number of rows/elements as x
#' @param ... additional arguments passed to methods
#' @return A `multiframe` object
#' @seealso [multidesign()]
#' @examples
#' # Create sample data
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble(condition = rep(letters[1:5], 4))
#'
#' # Create multiframe object
#' mf <- multiframe(X, Y)
#' @export
multiframe <- function(x, y, ...) UseMethod("multiframe")




#' Create a Multiblock Object
#'
#' Construct a new multiblock matrix consisting of a set of stacked submatrices sharing 
#' a row or column dimension. This structure is useful for analyzing data with multiple 
#' related blocks of measurements.
#'
#' @param x a list of matrices sharing either row or column dimensions
#' @param ... additional arguments passed to methods
#' @return A `multiblock` object
#' @seealso [is_cstacked()], [is_rstacked()]
#' @examples
#' # Create list of matrices with varying row dimensions
#' X <- lapply(c(10,20,30), function(i) matrix(rnorm(i*100), i, 100))
#' mb <- multiblock(X)
#' @export
multiblock <- function(x, ...) UseMethod("multiblock")


#' Test if Multiblock Object is Column Stacked
#'
#' @param x the multiblock object to test
#' @return logical indicating if the object is column stacked
#' @seealso [is_rstacked()], [multiblock()]
#' @export
is_cstacked <- function(x) UseMethod("is_cstacked")


#' Test if Multiblock Object is Row Stacked
#'
#' @param x the multiblock object to test
#' @return logical indicating if the object is row stacked
#' @seealso [is_cstacked()], [multiblock()]
#' @export
is_rstacked <- function(x) UseMethod("is_rstacked")


#' Summarize Data Over Grouping Variables
#'
#' Apply a multivariate columnwise summary function (e.g. `colMeans`) to subsets
#' of observations formed by one or more grouping variables.
#'
#' @param x the multivariate data (a matrix, list, or other data container)
#' @param sfun the columnwise summary function (e.g. `colMeans`)
#' @param ... the grouping variables
#' @return A summarized version of the input data
#' @examples
#' X <- matrix(rnorm(100*10), 100, 10)
#' groups <- rep(letters[1:5], each=20)
#' summarize_by(X, colMeans, groups)
#' @export
summarize_by <- function(x, sfun, ...) UseMethod("summarize_by")


#' Generate Cross-validation Folds
#'
#' Create a set of cross-validation objects over a blocking variable (`stratum`).
#' This function helps in creating stratified cross-validation splits while 
#' respecting the structure of the data.
#'
#' @param x the dataset to fold over
#' @param ... additional arguments passed to methods, such as number of folds or stratification variables
#' @return A list of fold indices for cross-validation
#' @examples
#' # Create example data with 100 observations and 2 groups
#' X <- matrix(rnorm(100*10), 100, 10)
#' groups <- rep(1:2, each=50)
#' folds <- fold_over(X, nfolds=5, stratum=groups)
#' @export
fold_over <- function(x, ...) UseMethod("fold_over")

#' Extract Data Matrix
#'
#' Get the data matrix (`x`) from a multidesign or related object.
#'
#' @param x the object containing the data matrix
#' @param ... additional arguments passed to methods
#' @return The data matrix component of the object
#' @seealso [design()]
#' @export
xdata <- function(x, ...) UseMethod("xdata")

#' Extract Design Information
#'
#' Get design variables associated with the rows of a multidesign or related object.
#'
#' @param x the object containing design information
#' @param ... additional arguments passed to methods
#' @return The design component (typically a data frame or matrix) of the object
#' @seealso [xdata()], [column_design()]
#' @export
design <- function(x, ...) UseMethod("design")

#' Extract Column Design Information
#'
#' Get design variables associated with the columns of a multidesign or related object.
#' This is particularly useful for datasets where variables (columns) have associated metadata.
#'
#' @param x the object containing column design information
#' @param ... additional arguments passed to methods
#' @return The column design component of the object
#' @seealso [design()]
#' @export
column_design <- function(x, ...) UseMethod("column_design")

#' Split Indices by Variables
#'
#' Extract the row indices of a table, split by one or more variables. This is useful
#' for creating grouped analyses or stratified sampling.
#'
#' @param x the object to split
#' @param ... the splitting variables
#' @return A list of integer vectors containing indices for each split
#' @examples
#' # Create example data
#' dat <- data.frame(x=1:100, 
#'                   group=rep(letters[1:4], each=25))
#' # Split by group
#' indices <- split_indices(dat, group)
#' @export
split_indices <- function(x, ...) UseMethod("split_indices")


#' Initialize a Transform Object
#'
#' @title Initialize a Transform Object
#' @description This is a re-export of the init_transform function from the multivarious package.
#' See \code{multivarious::init_transform} for full documentation.
#' @name init_transform
#' @importFrom multivarious init_transform
#' @export
multivarious::init_transform

#' Get Block Indices
#'
#' Retrieve the indices that define the boundaries of blocks in a multiblock or hyperdesign object.
#'
#' @param x The object to get block indices from
#' @param ... Additional arguments passed to methods
#' @return A matrix with start and end indices for each block
#' @export
block_indices <- function(x, ...) UseMethod("block_indices")

#' Get Block Indices
#'
#' @title Get Block Indices from a Multiblock Object
#' @description This is a re-export of the block_indices function from the multivarious package.
#' See \code{multivarious::block_indices} for full documentation.
#' @name block_indices
#' @importFrom multivarious block_indices
#' @export
multivarious::block_indices
