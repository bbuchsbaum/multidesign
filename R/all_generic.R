# Suppress R CMD check notes for NSE variables
utils::globalVariables(c("data", "temp"))

#' Create an Observation Object
#'
#' @description
#' Constructs a new lazy-evaluated observation object from various data sources.
#' An observation represents a single row, element, or vector from a data source
#' that can be accessed on demand.
#'
#' @param x The data source (matrix, list, vector, or other supported object)
#' @param i The index of the observation to extract
#'
#' @return A function that, when called, returns the specified observation from the data source
#'
#' @examples
#' # From a matrix
#' X <- matrix(1:20, 5, 4)
#' obs1 <- observation(X, 2)  # Create observation for row 2
#' obs1()  # Retrieve the observation (returns row 2 as a 1xn matrix)
#'
#' # From a list
#' X_list <- list(a=1:5, b=6:10, c=11:15)
#' obs2 <- observation(X_list, 3)  # Create observation for 3rd element
#' obs2()  # Retrieve the observation (returns the 3rd element: 11:15)
#'
#' @seealso [multiframe()], [obs_group()]
#' @export
observation <- function(x, i) UseMethod("observation")



#' Create a Multidesign Object
#'
#' @description
#' Constructs a new multivariate design object linking vector-valued observations with
#' design variables. A multidesign object maintains the relationship between experimental
#' data (observations) and metadata about experimental conditions (design variables).
#'
#' @details
#' A multidesign object consists of three main components:
#' * A data matrix where rows represent observations and columns represent variables
#' * A design data frame containing experimental factors and conditions for each observation
#' * Optional column metadata describing properties of each variable
#'
#' @param x The multivariate data (a matrix, a list, or other data container)
#' @param y A design matrix or data frame with same number of rows/elements as x
#' @param ... Additional arguments passed to methods, such as column_design
#'
#' @return A multidesign object with components:
#'   \item{x}{The input data matrix}
#'   \item{design}{A tibble containing design variables}
#'   \item{column_design}{A tibble containing column metadata (if provided)}
#'
#' @examples
#' # Create example data matrix
#' X <- matrix(rnorm(20*100), 20, 100)
#' 
#' # Create design information
#' Y <- tibble::tibble(condition=rep(letters[1:5], 4))
#' 
#' # Create multidesign object
#' mds <- multidesign(X, Y)
#' 
#' # Split by condition
#' sdes <- split(mds, condition)
#'
#' @seealso 
#'   \code{\link{reduce.multidesign}} for dimensionality reduction,
#'   \code{\link{split.multidesign}} for splitting by design variables,
#'   \code{\link{multiframe}} for an alternative implementation
#' @export
multidesign <- function(x, y, ...) UseMethod("multidesign")

#' Create a Hyperdesign Object
#'
#' @description
#' Constructs a new `hyperdesign` object that encapsulates a collection of `multidesign` objects.
#' Used to model multi-group, multi-block, or multi-view datasets, where each block/group/view
#' is associated with a matrix-variate response and an arbitrary design.
#'
#' @details
#' A hyperdesign object represents a collection of related multivariate datasets (multidesign instances)
#' that share common design variables. This structure is particularly useful for:
#' * Multiple subjects in an experiment
#' * Multiple sessions or runs
#' * Multiple data modalities (e.g., fMRI, EEG, behavioral)
#' * Multiple response measures
#'
#' @param x A list of multidesign objects. Each instance should represent a related block of data
#' @param block_names Optional character vector of names for each block
#'
#' @return A hyperdesign object with the following components:
#'   \item{blocks}{List of multidesign objects}
#'   \item{block_names}{Names of each block}
#'   \item{col_indices}{Matrix of column start/end indices for each block}
#'   \item{row_indices}{Matrix of row start/end indices for each block}
#'
#' @examples
#' # Create three multidesign objects (e.g., for three subjects)
#' d1 <- multidesign(
#'   matrix(rnorm(10*20), 10, 20),
#'   data.frame(y=1:10, subject=1, run=rep(1:5, 2))
#' )
#' d2 <- multidesign(
#'   matrix(rnorm(10*20), 10, 20),
#'   data.frame(y=1:10, subject=2, run=rep(1:5, 2))
#' )
#' d3 <- multidesign(
#'   matrix(rnorm(10*20), 10, 20),
#'   data.frame(y=1:10, subject=3, run=rep(1:5, 2))
#' )
#'
#' # Combine into a hyperdesign
#' hd <- hyperdesign(
#'   list(d1, d2, d3),
#'   block_names = c("subject1", "subject2", "subject3")
#' )
#'
#' @seealso 
#'   \code{\link{df_to_hyperdesign}} for creating hyperdesign objects from data frames,
#'   \code{\link{multidesign}} for the underlying multidesign structure,
#'   \code{\link{multiblock}} for another multi-block data structure
#' @export
hyperdesign <- function(x, block_names = NULL) UseMethod("hyperdesign")



#' Create a Multiframe Object
#'
#' @description
#' Constructs a new multivariate design object linking vector-valued observations with
#' arbitrary design variables. A multiframe combines experimental design metadata with 
#' lazy-evaluated observations, providing a flexible interface for data manipulation.
#'
#' @details
#' A multiframe object is similar to a multidesign but with enhanced support for data frame 
#' operations. Key features include:
#' * Lazy evaluation of observations (only loaded when needed)
#' * Integration with tidyverse functions for data manipulation
#' * Support for various data sources (matrices, lists, vectors)
#' * Methods for splitting, summarizing, and transforming data
#'
#' @param x The multivariate data (a matrix, list, or other data container)
#' @param y A design matrix or data frame with same number of rows/elements as x
#' @param ... Additional arguments passed to methods
#'
#' @return A multiframe object containing:
#'   \item{design}{A tibble with design variables and observation functions}
#'
#' @examples
#' # Create sample data
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble::tibble(condition = rep(letters[1:5], 4))
#'
#' # Create multiframe object
#' mf <- multiframe(X, Y)
#' 
#' # Access first observation
#' obs1 <- mf$design$.obs[[1]]()
#' 
#' # Split by condition
#' split_by_cond <- split(mf, condition)
#' 
#' # Summarize by condition
#' means_by_cond <- summarize_by(mf, condition)
#'
#' @seealso 
#'   \code{\link{multidesign}} for an alternative implementation,
#'   \code{\link{observation}} for the underlying observation structure,
#'   \code{\link{obs_group}} for creating groups of observations
#' @export
multiframe <- function(x, y, ...) UseMethod("multiframe")




#' Create a Multiblock Object
#'
#' @description
#' Constructs a new multiblock object consisting of a set of stacked submatrices sharing
#' a row or column dimension. This structure is useful for analyzing data with multiple
#' related blocks of measurements while preserving block structure information.
#'
#' @details
#' A multiblock object automatically determines whether matrices should be:
#' * Row-stacked: all matrices must have the same number of columns
#' * Column-stacked: all matrices must have the same number of rows
#'
#' The resulting object maintains block structure information while allowing operations
#' across the entire combined matrix.
#'
#' @param x A list of matrices sharing either row or column dimensions
#' @param ... Additional arguments passed to methods
#'
#' @return A multiblock object with the following attributes:
#'   \item{ind}{Matrix of start/end indices for each block}
#'   \item{orient}{Orientation of stacking ("cstacked" or "rstacked")}
#'
#' @examples
#' # Create list of matrices with varying row dimensions (column-stacked)
#' X1 <- matrix(rnorm(20*3), 20, 3)
#' X2 <- matrix(rnorm(20*5), 20, 5)
#' X3 <- matrix(rnorm(20*4), 20, 4)
#' mb_c <- multiblock(list(X1, X2, X3))
#' is_cstacked(mb_c)  # TRUE
#'
#' # Create list of matrices with varying column dimensions (row-stacked)
#' Y1 <- matrix(rnorm(10*5), 10, 5)
#' Y2 <- matrix(rnorm(15*5), 15, 5)
#' Y3 <- matrix(rnorm(20*5), 20, 5)
#' mb_r <- multiblock(list(Y1, Y2, Y3))
#' is_rstacked(mb_r)  # TRUE
#'
#' @seealso 
#'   \code{\link{is_cstacked}} for checking if a multiblock object is column-stacked,
#'   \code{\link{is_rstacked}} for checking if a multiblock object is row-stacked,
#'   \code{\link{block_indices}} for accessing block-specific indices
#' @export
multiblock <- function(x, ...) UseMethod("multiblock")


#' Test if Multiblock Object is Column Stacked
#'
#' @description
#' Checks if a multiblock object is column-stacked, meaning that all component matrices
#' share the same number of rows. In column-stacked multiblock objects, blocks are arranged
#' side by side, with each block potentially having a different number of columns.
#'
#' @param x The multiblock object to test
#' @return Logical value: TRUE if the object is column-stacked, FALSE otherwise
#'
#' @examples
#' # Create column-stacked multiblock (matrices share row dimension)
#' X1 <- matrix(rnorm(20*3), 20, 3)
#' X2 <- matrix(rnorm(20*5), 20, 5)
#' mb <- multiblock(list(X1, X2))
#' 
#' # Test stacking orientation
#' is_cstacked(mb)  # Returns TRUE
#' is_rstacked(mb)  # Returns FALSE
#'
#' @seealso 
#'   \code{\link{is_rstacked}} for testing if a multiblock object is row-stacked,
#'   \code{\link{multiblock}} for creating multiblock objects
#' @export
is_cstacked <- function(x) UseMethod("is_cstacked")


#' Test if Multiblock Object is Row Stacked
#'
#' @description
#' Checks if a multiblock object is row-stacked, meaning that all component matrices
#' share the same number of columns. In row-stacked multiblock objects, blocks are arranged
#' one above another, with each block potentially having a different number of rows.
#'
#' @param x The multiblock object to test
#' @return Logical value: TRUE if the object is row-stacked, FALSE otherwise
#'
#' @examples
#' # Create row-stacked multiblock (matrices share column dimension)
#' Y1 <- matrix(rnorm(10*5), 10, 5)
#' Y2 <- matrix(rnorm(15*5), 15, 5)
#' mb <- multiblock(list(Y1, Y2))
#' 
#' # Test stacking orientation
#' is_rstacked(mb)  # Returns TRUE
#' is_cstacked(mb)  # Returns FALSE
#'
#' @seealso 
#'   \code{\link{is_cstacked}} for testing if a multiblock object is column-stacked,
#'   \code{\link{multiblock}} for creating multiblock objects
#' @export
is_rstacked <- function(x) UseMethod("is_rstacked")


#' Summarize Data Over Grouping Variables
#'
#' @description
#' Applies a multivariate columnwise summary function (e.g., `colMeans`) to subsets
#' of observations formed by one or more grouping variables. This function provides
#' a flexible way to compute group-wise statistics across different data structures.
#'
#' @details
#' The function works with various data structures including:
#' * Multidesign objects: Summarizes the data matrix by design variables
#' * Multiframe objects: Summarizes lazy-evaluated observations by design variables
#' * Matrices: Summarizes columns by external grouping variables
#'
#' For each unique combination of grouping variables, the summary function is applied
#' to the corresponding subset of observations.
#'
#' @param x The multivariate data (a matrix, multidesign, multiframe, or other data container)
#' @param ... Method-specific arguments. For multidesign/multiframe: unquoted names of variables
#'   to group by, plus optional `sfun` (summary function, default `colMeans`) and
#'   `extract_data` (logical). For matrices: the grouping variables.
#'
#' @return For multidesign/multiframe objects: a new object of the same class with summarized data.
#'         For matrices: a matrix of summary statistics for each group.
#'
#' @examples
#' # With a multidesign object
#' X <- matrix(rnorm(100*10), 100, 10)
#' groups <- rep(letters[1:5], each=20)
#' mds <- multidesign(X, data.frame(group=groups))
#' mds_means <- summarize_by(mds, group)
#'
#' # With a multiframe object
#' mf <- multiframe(X, data.frame(group=groups))
#' mf_means <- summarize_by(mf, group)
#'
#' @seealso
#'   \code{\link{multidesign}} for creating multidesign objects,
#'   \code{\link{multiframe}} for creating multiframe objects,
#'   \code{\link{split}} for splitting data by grouping variables
#' @export
summarize_by <- function(x, ...) UseMethod("summarize_by")


#' Generate Cross-validation Folds
#'
#' @description
#' Creates cross-validation folds by splitting data based on specified design variables
#' or stratification factors. This function helps in creating stratified cross-validation
#' splits while respecting the structure of the data, particularly for complex experimental
#' designs.
#'
#' @details
#' The function creates folds by splitting the data based on unique combinations of the specified
#' variables. For each fold, one combination is held out as the assessment set (test set), 
#' while the rest form the analysis set (training set).
#'
#' Different behaviors are implemented for different object types:
#' * For multidesign objects: Creates folds based on design variables
#' * For hyperdesign objects: Creates folds across blocks, respecting block structure
#' * For matrices: Creates folds based on external stratification variables
#'
#' @param x The dataset to fold over (multidesign, hyperdesign, matrix, etc.)
#' @param ... Additional arguments passed to methods:
#'   * For multidesign/hyperdesign: unquoted names of variables to split on
#'   * For multidesign/hyperdesign/multiframe: `preserve_row_ids = TRUE` to carry
#'     original source row ids into fold outputs via a reserved `.orig_index` column
#'   * For matrices: nfolds (number of folds), stratum (stratification variable)
#'
#' @return A foldlist object containing:
#'   \item{analysis}{Training data for each fold}
#'   \item{assessment}{Test data for each fold}
#'   \item{held_out}{Information about which values were held out in each fold}
#'
#' @examples
#' # With a multidesign object
#' X <- matrix(rnorm(100*10), 100, 10)
#' groups <- rep(1:2, each=50)
#' mds <- multidesign(X, data.frame(group=groups, subject=rep(1:10, each=10)))
#' folds_by_group <- fold_over(mds, group)
#'
#' # With a hyperdesign object (multiple subjects)
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                  data.frame(condition=rep(c("A","B"), 5), run=1:10))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                  data.frame(condition=rep(c("A","B"), 5), run=1:10))
#' hd <- hyperdesign(list(d1, d2), block_names=c("subject1", "subject2"))
#' folds_by_condition <- fold_over(hd, condition)
#'
#' @seealso 
#'   \code{\link{multidesign}} for creating multidesign objects,
#'   \code{\link{hyperdesign}} for creating hyperdesign objects,
#'   \code{\link{split_indices}} for splitting indices without creating folds
#' @export
fold_over <- function(x, ...) UseMethod("fold_over")

#' Extract Data Matrix
#'
#' @description
#' Retrieves the data matrix component from various object types in the package.
#' This function provides a consistent interface for accessing the underlying data
#' regardless of the specific object structure.
#'
#' @details
#' The function behaves differently depending on the class of the input object:
#' * For multidesign objects: Returns the data matrix component
#' * For hyperdesign objects: Returns a list of data matrices, one for each block
#'   (or a single matrix if a specific block is requested)
#' * For multiframe objects: Returns the combined data from all observations
#'
#' @param x The object containing the data matrix
#' @param ... Additional arguments passed to methods:
#'   * For hyperdesign objects: 'block' parameter to specify which block's data to return
#'
#' @return The data matrix component of the object, or a list of matrices for hyperdesign objects
#'
#' @examples
#' # With a multidesign object
#' X <- matrix(rnorm(20*10), 20, 10)
#' Y <- data.frame(group = rep(letters[1:4], each=5))
#' mds <- multidesign(X, Y)
#' X_data <- xdata(mds)  # Returns the original matrix X
#'
#' # With a hyperdesign object
#' d1 <- multidesign(matrix(rnorm(10*5), 10, 5), data.frame(subject=rep(1,10)))
#' d2 <- multidesign(matrix(rnorm(10*5), 10, 5), data.frame(subject=rep(2,10)))
#' hd <- hyperdesign(list(d1, d2))
#' all_data <- xdata(hd)  # Returns list of matrices
#' block1_data <- xdata(hd, block=1)  # Returns just the first block's matrix
#'
#' @seealso 
#'   \code{\link{design}} for extracting design information,
#'   \code{\link{column_design}} for extracting column metadata
#' @export
xdata <- function(x, ...) UseMethod("xdata")

#' Extract Design Information
#'
#' @description
#' Retrieves the design information (experimental metadata) associated with observations
#' in various object types. This function provides a consistent interface for accessing
#' design variables regardless of the specific object structure.
#'
#' @details
#' The function behaves differently depending on the class of the input object:
#' * For multidesign objects: Returns the design data frame component
#' * For hyperdesign objects: Returns a list of design data frames, one for each block
#'   (or a single design if a specific block is requested)
#' * For multiframe objects: Returns the design tibble with observation functions
#'
#' Design information typically includes experimental factors, conditions, subject IDs,
#' or other metadata associated with each observation.
#'
#' @param x The object containing design information
#' @param ... Additional arguments passed to methods:
#'   * For hyperdesign objects: 'block' parameter to specify which block's design to return
#'
#' @return The design component of the object (typically a data frame or tibble),
#'         or a list of designs for hyperdesign objects
#'
#' @examples
#' # With a multidesign object
#' X <- matrix(rnorm(20*10), 20, 10)
#' Y <- data.frame(
#'   condition = rep(c("A", "B"), each=10),
#'   subject = rep(1:5, times=4)
#' )
#' mds <- multidesign(X, Y)
#' design_info <- design(mds)  # Returns the original data frame Y
#'
#' # With a hyperdesign object
#' d1 <- multidesign(matrix(rnorm(10*5), 10, 5), 
#'                  data.frame(condition=rep(c("A","B"), 5), subject=1))
#' d2 <- multidesign(matrix(rnorm(10*5), 10, 5), 
#'                  data.frame(condition=rep(c("A","B"), 5), subject=2))
#' hd <- hyperdesign(list(d1, d2))
#' all_designs <- design(hd)  # Returns list of design data frames
#' block1_design <- design(hd, block=1)  # Returns just the first block's design
#'
#' @seealso 
#'   \code{\link{xdata}} for extracting the data matrix,
#'   \code{\link{column_design}} for extracting column metadata
#' @export
design <- function(x, ...) UseMethod("design")

#' Extract Column Design Information
#'
#' @description
#' Retrieves metadata associated with the columns (variables) of a multidesign or related object.
#' This function is particularly useful for datasets where variables have associated metadata,
#' such as feature names, measurement types, or other variable-specific information.
#'
#' @details
#' The function behaves differently depending on the class of the input object:
#' * For multidesign objects: Returns the column_design data frame component
#' * For hyperdesign objects: Returns a list of column design data frames, one for each block
#'   (or a single column design if a specific block is requested)
#'
#' Column design information typically includes variable names, feature types, regions of interest,
#' or other metadata that describes the variables rather than the observations.
#'
#' @param x The object containing column design information
#' @param ... Additional arguments passed to methods:
#'   * For hyperdesign objects: 'block' parameter to specify which block's column design to return
#'
#' @return The column design component of the object (typically a data frame or tibble),
#'         or a list of column designs for hyperdesign objects
#'
#' @examples
#' # With a multidesign object including column metadata
#' X <- matrix(rnorm(20*10), 20, 10)
#' Y <- data.frame(condition = rep(c("A", "B"), each=10))
#' col_info <- data.frame(
#'   feature = paste0("var", 1:10),
#'   type = rep(c("continuous", "categorical"), 5)
#' )
#' mds <- multidesign(X, Y, col_info)
#' col_metadata <- column_design(mds)
#'
#' # With a hyperdesign object
#' d1 <- multidesign(matrix(rnorm(10*5), 10, 5), 
#'                  data.frame(condition=rep(c("A","B"), 5)),
#'                  data.frame(feature=paste0("var", 1:5)))
#' d2 <- multidesign(matrix(rnorm(10*5), 10, 5), 
#'                  data.frame(condition=rep(c("A","B"), 5)),
#'                  data.frame(feature=paste0("var", 1:5)))
#' hd <- hyperdesign(list(d1, d2))
#' all_col_designs <- column_design(hd)
#' block1_col_design <- column_design(hd, block=1)
#'
#' @seealso 
#'   \code{\link{design}} for extracting observation design information,
#'   \code{\link{xdata}} for extracting the data matrix
#' @export
column_design <- function(x, ...) UseMethod("column_design")

#' Split Indices by Variables
#'
#' @description
#' Extracts the row indices of a dataset, split by one or more variables. This function
#' is useful for creating grouped analyses, stratified sampling, or preparing data for
#' cross-validation without actually duplicating the data.
#'
#' @details
#' The function returns indices rather than the actual data, which is memory-efficient
#' for large datasets. It works with various object types:
#' * For data frames: Returns indices grouped by specified variables
#' * For multidesign objects: Returns indices grouped by design variables
#' * For multiframe objects: Returns indices grouped by design variables
#'
#' The resulting indices can be used to subset the original data or to create
#' training/testing splits for machine learning.
#'
#' @param x The object to split (data frame, multidesign, multiframe, etc.)
#' @param ... Unquoted names of variables to split by
#' @param collapse Logical; whether to collapse the resulting indices
#'
#' @return A tibble or list containing:
#'   \item{group variables}{The splitting variables and their values}
#'   \item{indices}{List column of integer vectors with row indices for each group}
#'
#' @examples
#' # With a multidesign object
#' mds <- multidesign(matrix(rnorm(100*10), 100, 10),
#'                   data.frame(group=rep(letters[1:4], each=25),
#'                              block=rep(1:5, times=20)))
#' mds_indices <- split_indices(mds, group)
#'
#' @seealso 
#'   \code{\link{fold_over}} for creating cross-validation folds,
#'   \code{\link{split}} for splitting the actual data
#' @export
split_indices <- function(x, ...) UseMethod("split_indices")


#' Select Variables Based on Column Design
#'
#' @description
#' Subsets the columns (variables) of an object based on conditions applied to the
#' column design metadata. This allows filtering variables by their properties
#' (e.g., region, type, hemisphere) using dplyr-style filter expressions.
#'
#' @param x The object to select variables from (multidesign, hyperdesign, etc.)
#' @param ... Filter expressions evaluated against the column design (using dplyr::filter semantics)
#'
#' @return An object of the same class with a subset of columns/variables
#'
#' @examples
#' X <- matrix(rnorm(20*10), 20, 10)
#' Y <- data.frame(condition = rep(c("A", "B"), each=10))
#' col_info <- data.frame(
#'   region = rep(c("frontal", "parietal"), 5),
#'   hemisphere = rep(c("left", "right"), each=5)
#' )
#' mds <- multidesign(X, Y, col_info)
#'
#' # Select only frontal variables
#' mds_frontal <- select_variables(mds, region == "frontal")
#'
#' @seealso
#'   \code{\link{column_design}} for extracting column metadata,
#'   \code{\link{multidesign}} for creating multidesign objects
#' @export
select_variables <- function(x, ...) UseMethod("select_variables")

#' Convert to a Multidesign Object
#'
#' @description
#' Converts an object (e.g., a hyperdesign) into a single multidesign object
#' by collapsing its structure.
#'
#' @param x The object to convert
#' @param ... Additional arguments passed to methods
#'
#' @return A multidesign object
#'
#' @examples
#' d1 <- multidesign(matrix(rnorm(10*5), 10, 5),
#'                   data.frame(condition = rep(c("A","B"), 5)))
#' d2 <- multidesign(matrix(rnorm(10*5), 10, 5),
#'                   data.frame(condition = rep(c("A","B"), 5)))
#' hd <- hyperdesign(list(d1, d2))
#'
#' # Collapse hyperdesign to a single multidesign
#' md <- as_multidesign(hd)
#'
#' @seealso
#'   \code{\link{hyperdesign}} for creating hyperdesign objects,
#'   \code{\link{multidesign}} for the multidesign class
#' @export
as_multidesign <- function(x, ...) UseMethod("as_multidesign")

#' Initialize Transformation
#'
#' @title Initialize Transformation
#' @description Re-exported from \pkg{multivarious}.
#' See \code{multivarious::init_transform} for details.
#'
#' @param x The object to transform
#' @param X The data or preprocessing specification
#' @param ... Additional arguments passed to methods
#' @return The transformed object (class-specific)
#' @name init_transform
#' @importFrom multivarious init_transform
#' @export
multivarious::init_transform

#' Get Block Indices
#'
#' @title Get Block Indices from a Multiblock Object
#' @description This is a re-export of the block_indices function from the multivarious package.
#' See \code{multivarious::block_indices} for full documentation.
#' @param x The object to get block indices from
#' @param ... Additional arguments passed to methods
#' @return An integer vector or list of integer vectors with indices for the
#'   requested block(s)
#'
#' @examples
#' mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
#' block_indices(mb, 1)
#' block_indices(mb, 2)
#'
#' @name block_indices
#' @importFrom multivarious block_indices
#' @export
multivarious::block_indices
