#' Create a Multiframe from a List
#'
#' @description
#' Creates a multiframe object from a list of observations and associated design information.
#' A multiframe combines experimental design metadata with lazy-evaluated observations.
#' Each element in the list must have consistent dimensions (same number of columns if matrices,
#' or same length if vectors).
#'
#' @param x A list containing observations (matrices, arrays, or vectors)
#' @param y A data frame containing design variables. Must have same number of rows as length of x
#' @param ... Additional arguments (currently unused)
#' @return A multiframe object containing a design tibble with observation functions
#'
#' @examples
#' # Create list of observations (matrices with same number of columns)
#' x <- list(
#'   matrix(1:12, 3, 4),  # 3x4 matrix
#'   matrix(13:24, 3, 4), # 3x4 matrix
#'   matrix(25:36, 3, 4)  # 3x4 matrix
#' )
#' 
#' # Create design information
#' y <- data.frame(
#'   condition = c("A", "B", "C"),
#'   block = 1:3
#' )
#' 
#' # Create multiframe
#' mf <- multiframe(x, y)
#' 
#' # Access first observation
#' obs1 <- mf$design$.obs[[1]]()
#' 
#' # View design information
#' print(mf$design)
#'
#' @family multiframe functions
#' @export
multiframe.list <- function(x, y, ...) {
  chk::chk_list(x)
  chk::chk_not_empty(x)
  
  # Check if all elements are matrices/arrays/vectors
  valid_types <- sapply(x, function(elem) {
    is.matrix(elem) || is.array(elem) || is.vector(elem)
  })
  if (!all(valid_types)) {
    stop("All elements in x must be matrices, arrays, or vectors")
  }
  
  # Get number of rows and columns for each element
  dims <- lapply(x, function(elem) {
    if (is.matrix(elem) || is.array(elem)) {
      c(nrow(elem), ncol(elem))
    } else if (is.vector(elem)) {
      if (length(elem) == 0) {
        c(0, 0)
      } else {
        c(1, length(elem))
      }
    }
  })
  
  # Check for empty elements
  empty_elems <- sapply(dims, function(d) any(d == 0))
  if (any(empty_elems)) {
    stop("All elements must be non-empty")
  }
  
  # Check if all elements have the same number of columns
  ncols <- sapply(dims, `[`, 2)
  if (length(unique(ncols)) > 1) {
    stop("All elements must have the same number of columns/length. Found dimensions: ", 
         paste(ncols, collapse=", "))
  }
  
  # Check if all elements have the same number of rows
  nrows <- sapply(dims, `[`, 1)
  if (length(unique(nrows)) > 1) {
    stop("All elements must have the same number of rows. Found dimensions: ", 
         paste(nrows, collapse=", "))
  }

  y <- tibble::as_tibble(y)
  chk::chk_equal(length(x), nrow(y))
  
  # Create indices first
  indices <- seq_len(nrow(y))
  
  # Create observation functions with captured indices
  obs_fns <- lapply(indices, function(i) {
    observation.list(x, i)
  })
  
  # Add to design tibble
  des <- y
  des$.index <- indices
  des$.obs <- obs_fns
  
  structure(list(design=des), class="multiframe")
}

#' Create a Multiframe from a Matrix
#'
#' @description
#' Creates a multiframe object from a matrix of observations and associated design information.
#' Each row of the matrix becomes a lazy-evaluated observation.
#'
#' @param x A matrix where rows are observations and columns are variables
#' @param y A data frame containing design variables. Must have same number of rows as x
#' @param ... Additional arguments (currently unused)
#' @return A multiframe object containing a design tibble with observation functions
#'
#' @examples
#' # Create matrix of observations (10 observations x 4 variables)
#' x <- matrix(1:40, 10, 4)
#' colnames(x) <- paste0("var", 1:4)  # Optional: name the variables
#' 
#' # Create design information with two factors
#' y <- data.frame(
#'   condition = rep(c("A", "B"), each=5),
#'   subject = rep(1:5, times=2)
#' )
#' 
#' # Create multiframe
#' mf <- multiframe(x, y)  # Using the generic method
#' 
#' # Access first observation
#' first_obs <- mf$design$.obs[[1]]()
#' print(first_obs)  # Shows the first row of x as a matrix
#' 
#' # View design structure with observation functions
#' head(mf$design)
#'
#' @family multiframe functions
#' @importFrom dplyr mutate rowwise n
#' @export
multiframe.matrix <- function(x, y, ...) {
  chk::chk_matrix(x)
  y <- tibble::as_tibble(y)
  chk::chk_equal(nrow(x), nrow(y))
  
  # Create indices first
  indices <- seq_len(nrow(y))
  
  # Create observation functions with captured indices
  obs_fns <- lapply(indices, function(i) {
    observation.matrix(x, i)
  })
  
  # Add to design tibble
  des <- y
  des$.index <- indices
  des$.obs <- obs_fns
  
  structure(list(design=des), class="multiframe")
}

#' Create an Observation Group
#'
#' @title Create a Group of Observations from Matrix or List Data
#'
#' @description
#' Creates an observation group object that provides a unified interface for accessing
#' individual observations from different data structures (matrices, lists, or deflists).
#' Each observation is wrapped in a function that provides lazy evaluation.
#'
#' @param X A matrix, list, or deflist containing the data
#' @param fun Optional function to apply to each observation
#' @param ind Optional vector of indices for the observations. If NULL, uses sequential indices
#'
#' @return An observation_set object containing lazy-evaluated observations
#'
#' @examples
#' # From a matrix
#' X <- matrix(rnorm(100), 20, 5)
#' obs <- obs_group(X)
#' 
#' # From a list
#' X_list <- list(a=1:5, b=6:10, c=11:15)
#' obs_list <- obs_group(X_list)
#'
#' @family multiframe functions
#' @export
obs_group <- function(X, fun=NULL, ind=NULL) {
  chk::chk_not_empty(X)
  chk::chkor_vld(chk::vld_s3_class(X, "matrix"), 
                 chk::vld_s3_class(X, "list"), 
                 chk::vld_s3_class(X, "deflist"))

  ret <- if (inherits(X, "list") || inherits(X, "deflist")) {
    if (is.null(ind)) {
      ind <- 1:length(X)
    } else {
      chk::chk_equal(length(X), length(ind))
    }
    lapply(1:length(X), function(i) {
      observation.list(X, ind[i])
    })
  } else {
    if (is.null(ind)) {
      ind <- 1:nrow(X)
    } else {
      chk::chk_equal(length(X), length(ind))
    }
    lapply(1:nrow(X), function(i) {
      observation.matrix(X, ind[i])
    })
  }

  structure(ret, class="observation_set")
}

#' Extract Single Observation from Observation Set
#'
#' @description
#' Extracts a single observation from an observation_set using double bracket indexing.
#' The result is automatically evaluated.
#'
#' @param x An observation_set object
#' @param i Index of the observation to extract
#' @return The evaluated observation
#'
#' @family multiframe functions
#' @export
`[[.observation_set` <- function(x,i) {
  z <- NextMethod()
  z()
}

#' Extract Multiple Observations from Observation Set
#'
#' @description
#' Extracts multiple observations from an observation_set using single bracket indexing.
#' All selected observations are automatically evaluated.
#'
#' @param x An observation_set object
#' @param i Indices of observations to extract
#' @return List of evaluated observations
#'
#' @family multiframe functions
#' @export
`[.observation_set` <- function(x, i) {
  z <- NextMethod()
  lapply(z, function(zi) zi())
}

#' Create Observation from Deflist
#'
#' @description
#' Creates a lazy-evaluated observation from a deflist object.
#'
#' @param x A deflist object
#' @param i Index of the observation to create
#' @return An observation object (function) that returns the data when called
#'
#' @family multiframe functions
#' @export
observation.deflist <- function(x, i) {
  function() {
    x[[i]]
  }
}

#' Create Observation from Vector
#'
#' @description
#' Creates a lazy-evaluated observation from a vector.
#'
#' @param x A vector
#' @param i Index (ignored for vectors, as they represent a single observation)
#' @return An observation object (function) that returns the vector when called
#'
#' @family multiframe functions
#' @export
observation.vector <- function(x,i) {
  chk::chk_scalar(i)
  f <- function() {
    x
  }

  structure(f, i=i, class="observation")
}

#' Create Observation from Matrix
#'
#' @description
#' Creates a lazy-evaluated observation from a row of a matrix.
#'
#' @param x A matrix
#' @param i Row index of the observation to create
#' @return An observation object (function) that returns the matrix row when called
#'
#' @family multiframe functions
#' @export
observation.matrix <- function(x, i) {
  function() {
    matrix(x[i,], nrow=1)  # Ensure output is a 1xn matrix
  }
}

#' Create Observation from List
#'
#' @description
#' Creates a lazy-evaluated observation from a list element.
#'
#' @param x A list
#' @param i Index of the list element to create as an observation
#' @return An observation object (function) that returns the list element when called
#'
#' @family multiframe functions
#' @export
observation.list <- function(x, i) {
  function() {
    x[[i]]
  }
}

#' Split a Multiframe by Design Variables
#'
#' @description
#' Splits a multiframe into groups based on one or more design variables.
#' Optionally collapses the observations within each group into a single matrix.
#'
#' @param x A multiframe object
#' @param ... Unquoted names of variables to split by
#' @param collapse Logical; if TRUE, combines observations in each group into a matrix
#' @return A nested tibble containing the split data
#'
#' @examples
#' x <- matrix(1:40, 10, 4)
#' y <- data.frame(
#'   condition = rep(c("A", "B"), each=5),
#'   subject = rep(1:5, times=2)
#' )
#' mf <- multiframe(x, y)
#' 
#' # Split by condition
#' split_by_cond <- split(mf, condition)
#' 
#' # Split by condition and subject, collapsing observations
#' split_by_both <- split(mf, condition, subject, collapse=TRUE)
#'
#' @family multiframe functions
#' @export
split.multiframe <- function(x, ..., collapse=FALSE) {
  nest.by <- rlang::quos(...)
  ret <- x$design %>% nest_by(!!!nest.by)

  if (collapse) {
    ret <- ret %>% rowwise() %>% mutate(data={
      list(do.call(rbind, lapply(data$.obs, function(o) o())))
    })
  }

  ret
}

#' Summarize a Multiframe by Design Variables
#'
#' @description
#' Computes summaries of observations grouped by design variables.
#'
#' @param x A multiframe object
#' @param sfun Summary function to apply (default is colMeans)
#' @param extract_data Logical; if TRUE, returns just the summarized data
#' @param ... Unquoted names of variables to group by
#' @return A tibble containing summarized data for each group
#'
#' @examples
#' x <- matrix(1:40, 10, 4)
#' y <- data.frame(
#'   condition = rep(c("A", "B"), each=5),
#'   subject = rep(1:5, times=2)
#' )
#' mf <- multiframe(x, y)
#' 
#' # Get means by condition
#' means_by_cond <- summarize_by(mf, condition)
#'
#' @family multiframe functions
#' @export
summarize_by.multiframe <- function(x, ..., sfun=colMeans, extract_data=FALSE) {
  nest.by <- rlang::enquos(...)
  ret <- x$design %>% dplyr::group_by(!!!nest.by) %>% tidyr::nest()

  ret <- ret %>% dplyr::rowwise() %>% dplyr::mutate(data = {
    list(sfun(do.call(rbind, lapply(data$.obs, function(o) o()))))
  })

  if (extract_data) {
    ret <- do.call(rbind, ret$data %>% purrr::map(~ .x))
  }

  ret
}

#' Print Method for Observation Set Objects
#'
#' @description
#' Displays a formatted summary of an observation set object.
#'
#' @param x An observation_set object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#'
#' @method print observation_set
#' @export
print.observation_set <- function(x, ...) {
  n_obs <- length(x)
  
  cat(crayon::bold(crayon::blue("\n=== Observation Set ===\n")))
  cat("\n", crayon::green(paste0("Number of Observations: ", n_obs)), "\n")
  
  # Sample a few observations if there are many
  if (n_obs > 0) {
    cat(crayon::bold("\nSample Observations:"), "\n")
    sample_size <- min(3, n_obs)
    sample_indices <- sort(sample(n_obs, sample_size))
    
    for (i in sample_indices) {
      obs <- x[[i]]
      if (is.matrix(obs)) {
        dims <- dim(obs)
        cat("  ", crayon::white("-"), " Observation ", i, ": ", 
            crayon::green(paste0(dims[1], " x ", dims[2], " matrix")), "\n", sep="")
      } else {
        cat("  ", crayon::white("-"), " Observation ", i, ": ", 
            crayon::green(paste0(typeof(obs), " [length: ", length(obs), "]")), "\n", sep="")
      }
    }
    if (n_obs > sample_size) {
      cat("  ", crayon::white("-"), " ", crayon::yellow("... and", n_obs - sample_size, "more observations"), "\n")
    }
  }
  
  cat(crayon::bold(crayon::blue("\n=====================\n")))
  invisible(x)
}

#' Print Method for Multiframe Objects
#'
#' @description
#' Displays a formatted summary of a multiframe object.
#'
#' @param x A multiframe object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#'
#' @method print multiframe
#' @export
print.multiframe <- function(x, ...) {
  cat(crayon::bold(crayon::blue("\n=== Multiframe Object ===\n")))
  
  # Basic information
  n_obs <- nrow(x$design)
  design_vars <- names(x$design)[!names(x$design) %in% c(".index", ".obs")]
  
  cat("\n", crayon::green(paste0("Number of Observations: ", n_obs)), "\n")
  
  # Design variables
  cat(crayon::bold("\nDesign Variables:"), "\n")
  for (var in design_vars) {
    unique_vals <- unique(x$design[[var]])
    n_unique <- length(unique_vals)
    sample_vals <- if (n_unique > 5) {
      paste0(paste(head(unique_vals, 3), collapse=", "), 
             "...", 
             paste(tail(unique_vals, 2), collapse=", "))
    } else {
      paste(unique_vals, collapse=", ")
    }
    cat("  ", crayon::white("-"), " ", var, ": ", 
        crayon::green(n_unique), " levels (", sample_vals, ")\n", sep="")
  }
  
  # Sample observations
  if (n_obs > 0) {
    cat(crayon::bold("\nSample Observations:"), "\n")
    sample_size <- min(3, n_obs)
    sample_indices <- sort(sample(n_obs, sample_size))
    
    for (i in sample_indices) {
      obs <- x$design$.obs[[i]]()
      if (is.matrix(obs)) {
        dims <- dim(obs)
        cat("  ", crayon::white("-"), " Observation ", i, ": ", 
            crayon::green(paste0(dims[1], " x ", dims[2], " matrix")), "\n", sep="")
      } else {
        cat("  ", crayon::white("-"), " Observation ", i, ": ", 
            crayon::green(paste0(typeof(obs), " [length: ", length(obs), "]")), "\n", sep="")
      }
    }
    if (n_obs > sample_size) {
      cat("  ", crayon::white("-"), " ", crayon::yellow("... and", n_obs - sample_size, "more observations"), "\n")
    }
  }
  
  cat(crayon::bold(crayon::blue("\n===================\n")))
  invisible(x)
}
