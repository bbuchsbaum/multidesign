#' @keywords internal
block_index_mat <- function(x, byrow=FALSE) {
  lens <- if (byrow) {
    sapply(x, function(z) nrow(z$x))
  } else {
    sapply(x, function(z) ncol(z$x))
  }

  csum <- cumsum(lens)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}



#' Convert a Data Frame to a Hyperdesign Object
#'
#' This function creates a hyperdesign object from a data frame by splitting it according to a specified variable.
#' It's particularly useful when you have a wide-format data frame that needs to be converted into multiple
#' related multidesign objects, such as when dealing with multiple subjects or sessions in an experiment.
#'
#' @importFrom dplyr select everything
#' @importFrom dplyr arrange group_by nest_by
#' @importFrom purrr map
#'
#' @param data A data frame or tibble containing both design variables and response variables
#' @param design_vars Character vector specifying the names of design variables (e.g., conditions, factors)
#' @param x_vars Character vector specifying the names of response variables to extract
#' @param split_var Character string naming the variable to split the data on (e.g., "subject" or "session")
#'
#' @return A hyperdesign object containing multiple multidesign objects, one for each unique value in split_var
#' @seealso
#'   \code{\link{hyperdesign}} for creating hyperdesign objects directly from multidesign objects,
#'   \code{\link{multidesign}} for the underlying multidesign structure,
#'   \code{\link{design.hyperdesign}} for extracting design information
#'
#' @family hyperdesign functions
#' @export
#'
#' @examples
#' # Create a sample tibble with multiple subjects
#' sample_tibble <- tibble(
#'   felab = rep(1:2, each = 3),
#'   attention = rep(c("DA", "FA"), times = 3),
#'   basis = rep(c("basis01", "basis02", "basis03"), times = 2),
#'   subject = rep(1001:1002, each = 3),
#'   `1` = rnorm(6),  # response variable 1
#'   `2` = rnorm(6),  # response variable 2
#'   `3` = rnorm(6)   # response variable 3
#' )
#'
#' # Convert to hyperdesign, splitting by subject
#' hd <- df_to_hyperdesign(
#'   data = sample_tibble,
#'   design_vars = c("felab", "attention", "basis"),
#'   x_vars = as.character(1:3),
#'   split_var = "subject"
#' )
df_to_hyperdesign <- function(data, design_vars, x_vars, split_var) {
  # Arrange the data by the splitting variable
  data <- data %>% arrange(!!rlang::sym(split_var))

  # Create the nested dataframes
  nested_data <- data %>%
    group_by(!!rlang::sym(split_var)) %>%
    nest_by()

  # Convert each nested dataframe into a multidesign object
  multidesigns <- map(nested_data$data, function(df) {
    # Select X and design variables
    X <- select(df, all_of(x_vars))
    design <- select(df, all_of(design_vars))

    # Create multidesign object
    multidesign(as.matrix(X), design)
  })

  # Convert the list of multidesign objects into a hyperdesign
  hyperdesign(multidesigns, as.character(nested_data[[1]]))
}



#' Construct a Hyperdesign Object
#'
#' Creates a hyperdesign object, which represents a collection of related multivariate datasets
#' (multidesign instances) that share common design variables. This class is particularly useful
#' for modeling multi-block data, where you want to analyze multiple related matrices, such as:
#' * Multiple subjects in an experiment
#' * Multiple sessions or runs
#' * Multiple data modalities (e.g., fMRI, EEG, behavioral)
#' * Multiple response measures
#'
#' @param x A list of multidesign instances. Each instance should represent a related block of data
#' @param block_names Optional character vector of names for each block. If NULL, blocks will be
#'   automatically named as "block_1", "block_2", etc.
#'
#' @return A hyperdesign object with the following components:
#'   \item{blocks}{List of multidesign objects}
#'   \item{block_names}{Names of each block}
#'   \item{col_indices}{Matrix of column start/end indices for each block}
#'   \item{row_indices}{Matrix of row start/end indices for each block}
#'
#' @seealso
#'   \code{\link{df_to_hyperdesign}} for creating hyperdesign objects from data frames,
#'   \code{\link{multidesign}} for the underlying multidesign structure,
#'   \code{\link{design.hyperdesign}} for extracting design information
#'
#' @family hyperdesign functions
#' @export
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
hyperdesign <- function(x, block_names=NULL) {
  chk::chk_true(all(sapply(x, function(d) inherits(d, "multidesign"))))

  bind_col <- block_index_mat(x, byrow=FALSE)
  bind_row <- block_index_mat(x, byrow=TRUE)

  if (!is.null(block_names)) {
    chk::chk_true(length(block_names) == length(x))
    names(x) <- block_names
  } else if (is.null(names(x))) {
    block_names <- paste0("block_", 1:length(x))
    names(x) <- block_names
  }

  hdes <- lapply(1:length(x), function(i) {
    tibble(block=i, block_name=block_names[i],
           nr=nrow(x[[i]]$x), nxvar=ncol(x[[i]]$x), nyvar=ncol(x[[i]]$design),
           row_start=bind_row[i,1], row_end=bind_row[i,2],
           col_start=bind_col[i,1], col_end=bind_col[i,2])
  }) %>% bind_rows()

  cvars <- lapply(x, function(z) names(z$design))
  cvars <- Reduce(intersect, cvars)

  structure(x,
            hdes=hdes,
            common_vars=cvars,
            class="hyperdesign")
}

#' @export
#' @param i the block number
#' @param byrow if true, return row-oriented indices
#' @rdname block_indices
block_indices.hyperdesign <- function(x, i, byrow=FALSE) {
  hd <- attr(x, "hdes")

  if (missing(i)) {
    lapply(1:nrow(hd), function(j) {
      if (byrow) {
        seq(hd$row_start[j], hd$row_end[j])
      } else {
        seq(hd$col_start[j], hd$col_end[j])
      }
    })
  } else {
    chk::chk_range(i, c(1, nrow(hd)))
    if (byrow) {
      seq(hd$row_start[i], hd$row_end[i])
    } else {
      seq(hd$col_start[i], hd$col_end[i])
    }
  }
}

#' Initialize Transformation for Hyperdesign
#'
#' Method to initialize transformations (e.g., scaling, centering) for hyperdesign objects.
#' Each block in the hyperdesign gets its own transformation object.
#'
#' @param x A hyperdesign object
#' @param preproc A preprocessing specification (e.g., from recipes package)
#' @return A list of initialized preprocessing objects, one for each block
#' @family hyperdesign functions
#' @method init_transform hyperdesign
#' @export

init_transform.hyperdesign <- function(x, preproc) {
  ## pre-processors
  proclist <- lapply(seq_along(x), function(i) {
    multivarious:::fresh(preproc) %>% prep()
  })

  names(proclist) <- names(x)

  ## subject-split and pre-processed data
  out <- lapply(seq_along(proclist), function(i) {
    p <- proclist[[i]]
    Xi <- x[[i]]$x
    Xout <- multivarious::init_transform(p, Xi)
    multidesign(Xout, x[[i]]$design)
  })

  des <- hyperdesign(out)
  attr(des, "preproc") <- proclist
  des
}


#' Create Cross-Validation Folds from a Hyperdesign Object
#'
#' Creates cross-validation folds by splitting the data based on specified design variables.
#' Each fold consists of a training set (analysis) and a test set (assessment).
#'
#' @param x A hyperdesign object
#' @param ... Unquoted names of variables to split on (e.g., condition, subject, run)
#' @param inclusion_condition Optional list specifying values to include in assessment sets
#' @param exclusion_condition Optional list specifying values to exclude from assessment sets
#'
#' @details
#' The function creates folds by splitting the data based on unique combinations of the specified
#' variables. For each fold, one combination is held out as the assessment set, while the rest
#' form the analysis set.
#'
#' Important considerations:
#' * If a splitting variable is confounded with blocks (e.g., each subject is in a separate block),
#'   the function will fail as there would be no training data available for that block.
#' * Numeric variables (like run numbers) are handled by converting them to factors for splitting.
#' * The function preserves the design structure within each fold.
#'
#' @return A foldlist object containing the cross-validation folds
#' @export
fold_over.hyperdesign <- function(x, ..., inclusion_condition = list(), exclusion_condition = list()) {
  # Get the splitting variables properly
  dots <- rlang::enquos(...)
  split_vars <- sapply(dots, rlang::quo_name)

  # If no splitting variables provided, create leave-one-block-out folds
  if (length(split_vars) == 0) {
    # Create a fold for each block
    extract <- function(i) {
      # Block i is the test set
      test_design <- x[[i]]$design
      test_indices <- seq_len(nrow(test_design))

      # Create test dataset
      testdat <- multidesign(x[[i]]$x[test_indices, , drop = FALSE],
                            test_design[test_indices, , drop = FALSE],
                            x[[i]]$column_design)

      # Create training dataset from all other blocks
      traindat <- hyperdesign(x[-i])

      # Get held_out values
      held_out_values <- list(block = i)

      list(analysis = traindat,
           assessment = testdat,
           held_out = held_out_values)
    }

    # Create foldframe for block-wise folding
    foldframe <- tibble(
      .block = seq_along(x),
      indices = lapply(seq_along(x), function(i) seq_len(nrow(x[[i]]$design))),
      .splitvar = paste0("block_", seq_along(x)),
      .fold = seq_along(x)
    )
  } else {
    # Early detection of confounded variables
    for (var in split_vars) {
      # Get unique values for each block
      block_values <- lapply(seq_along(x), function(i) {
        vals <- x[[i]]$design[[var]]
        if (is.null(vals)) {
          stop("Variable '", var, "' not found in design of block ", i)
        }
        unique(vals)
      })

      # For each block, check if it has only one unique value
      for (i in seq_along(block_values)) {
        curr_vals <- block_values[[i]]
        if (length(curr_vals) == 1) {
          # Get all values from other blocks
          other_blocks_vals <- unique(unlist(block_values[-i]))

          # Check if this value appears in any other block
          if (!any(curr_vals %in% other_blocks_vals)) {
            stop(
              "Variable '", var, "' is confounded with blocks:\n",
              "  Block ", i, " contains only value '", curr_vals, "' which doesn't appear in other blocks.\n",
              "  This would result in no training data being available for this block.\n",
              "  Consider using a different splitting variable or restructuring your data."
            )
          }
        }
      }
    }

    # Validate inclusion_condition if it is non-empty
    if (length(inclusion_condition) > 0) {
      validate_inclusion_condition <- function(x, inclusion_condition) {
        all_design <- do.call(rbind, lapply(x, function(block) block$design))
        for (factor_name in names(inclusion_condition)) {
          if (!factor_name %in% colnames(all_design)) {
            stop(paste("Factor", factor_name, "not found in the design."))
          }
          if (!any(inclusion_condition[[factor_name]] %in% unique(all_design[[factor_name]]))) {
            stop(paste("Level", inclusion_condition[[factor_name]], "not found for factor", factor_name))
          }
        }
      }
      validate_inclusion_condition(x, inclusion_condition)
    }

    # Validate exclusion_condition if it is non-empty
    if (length(exclusion_condition) > 0) {
      validate_exclusion_condition <- function(x, exclusion_condition) {
        all_design <- do.call(rbind, lapply(x, function(block) block$design))
        for (factor_name in names(exclusion_condition)) {
          if (!factor_name %in% colnames(all_design)) {
            stop(paste("Factor", factor_name, "not found in the design."))
          }
          if (!any(exclusion_condition[[factor_name]] %in% unique(all_design[[factor_name]]))) {
            stop(paste("Level", exclusion_condition[[factor_name]], "not found for factor", factor_name))
          }
        }
      }
      validate_exclusion_condition(x, exclusion_condition)
    }

    splits <- lapply(seq_along(x), function(i) {
      d <- x[[i]]
      split_indices(d, ...) %>% mutate(.block = i)
    })

    foldframe <- splits %>% bind_rows() %>% mutate(.fold = 1:n())
  }

  condition_met <- function(i) {
    if (i > nrow(foldframe)) return(FALSE)

    block <- foldframe[[".block"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    # Ensure ind is numeric and valid
    if (length(ind) == 0) return(FALSE)

    # Convert indices to numeric if they're stored as character/factor
    ind <- as.integer(as.character(ind))
    if (any(is.na(ind))) return(FALSE)

    # Get the test design data
    test_design <- x[[block]]$design[ind, , drop = FALSE]
    if (nrow(test_design) == 0) return(FALSE)

    # Apply inclusion_condition
    if (length(inclusion_condition) > 0) {
      for (factor_name in names(inclusion_condition)) {
        test_design <- test_design[test_design[[factor_name]] %in% inclusion_condition[[factor_name]], , drop = FALSE]
      }
    }

    # Apply exclusion_condition
    if (length(exclusion_condition) > 0) {
      for (factor_name in names(exclusion_condition)) {
        test_design <- test_design[!test_design[[factor_name]] %in% exclusion_condition[[factor_name]], , drop = FALSE]
      }
    }

    # Check if there are any remaining rows after filtering
    nrow(test_design) > 0
  }

  indices <- which(sapply(seq_len(nrow(foldframe)), condition_met))

  if (length(indices) == 0) {
    stop("No valid folds created. Check inclusion and exclusion conditions.")
  }

  extract <- function(i) {
    block <- foldframe[[".block"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    # Convert indices to numeric if they're stored as character/factor
    ind <- as.integer(as.character(ind))

    test_design <- x[[block]]$design[ind, , drop = FALSE]

    # Apply inclusion_condition
    if (length(inclusion_condition) > 0) {
      for (factor_name in names(inclusion_condition)) {
        test_design <- test_design[test_design[[factor_name]] %in% inclusion_condition[[factor_name]], , drop = FALSE]
      }
    }

    # Apply exclusion_condition
    if (length(exclusion_condition) > 0) {
      for (factor_name in names(exclusion_condition)) {
        test_design <- test_design[!test_design[[factor_name]] %in% exclusion_condition[[factor_name]], , drop = FALSE]
      }
    }

    test_indices <- as.integer(rownames(test_design))
    testdat <- multidesign(x[[block]]$x[test_indices, , drop = FALSE],
                          test_design,
                          x[[block]]$column_design)

    traindat <- hyperdesign(lapply(seq_along(x), function(j) {
      if (j == block) {
        if (length(test_indices) == nrow(x[[j]]$x)) {
          stop("Block ", j, " contains only data for the assessment set.\n",
               "This typically happens when a splitting variable is confounded with blocks.\n",
               "Consider using a different splitting variable or restructuring your data.")
        }
        multidesign(x[[block]]$x[-test_indices, , drop = FALSE],
                   x[[block]]$design[-test_indices, , drop = FALSE],
                   x[[block]]$column_design)
      } else {
        x[[j]]
      }
    }))

    # Get held_out values from the filtered test_design
    held_out_values <- test_design %>%
      distinct() %>%
      slice(1) %>%
      as.list()

    list(analysis = traindat,
         assessment = testdat,
         held_out = held_out_values)
  }

  ret <- deflist(function(i) extract(indices[i]), len = length(indices))

  names(ret) <- paste0("fold_", seq_along(ret))
  class(ret) <- c("foldlist", class(ret))
  attr(ret, "foldframe") <- foldframe
  ret
}


#' Extract Design Information from Hyperdesign
#'
#' Retrieves design information from a hyperdesign object, either for all blocks
#' or for a specific block.
#'
#' @param x A hyperdesign object
#' @param block Optional numeric index specifying which block's design to return
#' @return If block is specified, returns design for that block; otherwise returns
#'         a list of designs for all blocks
#' @family hyperdesign functions
#' @method design hyperdesign
#' @export
#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(y=1:10, condition=rep(c("A","B"), 5)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(y=1:10, condition=rep(c("A","B"), 5)))
#' hd <- hyperdesign(list(d1, d2))
#'
#' # Get all designs
#' all_designs <- design(hd)
#'
#' # Get design for block 1
#' block1_design <- design(hd, block=1)
design.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, design)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    design(x[[block]])
  }
}



#' Extract Data Matrix from Hyperdesign
#'
#' Get the data matrix component from a hyperdesign object, either for all blocks
#' or for a specific block.
#'
#' @param x A hyperdesign object
#' @param block Optional numeric index specifying which block's data to return
#' @return If block is specified, returns data matrix for that block; otherwise returns
#'         a list of data matrices for all blocks
#' @family hyperdesign functions
#' @seealso [design.hyperdesign()], [column_design.hyperdesign()]
#' @method xdata hyperdesign
#' @export
#' @examples
#' # Create example hyperdesign
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(y=1:10, condition=rep(c("A","B"), 5)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(y=1:10, condition=rep(c("A","B"), 5)))
#' hd <- hyperdesign(list(d1, d2))
#'
#' # Get data from all blocks
#' all_data <- xdata(hd)
#'
#' # Get data from block 1
#' block1_data <- xdata(hd, block=1)
xdata.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, xdata)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    xdata(x[[block]])
  }
}

#' Get Column Design Information for Hyperdesign
#'
#' Extract column design information from a hyperdesign object.
#'
#' @param x A hyperdesign object
#' @param block Optional block index or name to get design for specific block
#' @param ... Additional arguments passed to methods
#' @return A list of column design information for each block
#' @method column_design hyperdesign
#' @export
column_design.hyperdesign <- function(x, block, ...) {
  if (missing(block)) {
    lapply(x, column_design)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    column_design(x[[block]])
  }
}

#' Subset a Hyperdesign Object
#'
#' Create a new hyperdesign object containing only the specified blocks.
#'
#' @param x A hyperdesign object
#' @param subset Logical vector indicating which blocks to keep
#' @param ... Additional arguments (not used)
#' @return A new hyperdesign object containing only the selected blocks
#' @family hyperdesign functions
#' @method subset hyperdesign
#' @export
#' @examples
#' # Create example hyperdesign
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(subject=1, condition=rep(c("A","B"), 5)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(subject=2, condition=rep(c("A","B"), 5)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20),
#'                   data.frame(subject=3, condition=rep(c("A","B"), 5)))
#' hd <- hyperdesign(list(d1, d2, d3))
#'
#' # Keep only blocks 1 and 3
#' subset_hd <- subset(hd, c(TRUE, FALSE, TRUE))
subset.hyperdesign <- function(x, fexpr) {
  out <- lapply(x, function(d) {
    subset(d, !!rlang::enquo(fexpr))
  })

  rem <- unlist(purrr::map(out, is.null))
  if (sum(rem) == length(x)) {
    stop("subset expression does not match any rows in hyperdesign `x`")
  }

  onam <- names(x)[!rem]
  hyperdesign(out[!rem], onam)
}

#' Print Method for Hyperdesign Objects
#'
#' Displays a detailed summary of a hyperdesign object, including information about each block's
#' structure, design variables, and dimensions.
#'
#' @param x A hyperdesign object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#' @method print hyperdesign
#' @export
print.hyperdesign <- function(x, ...) {
  # Header
  cat(bold(blue("\n═══ Hyperdesign Object ═══\n")))

  # Number of blocks
  cat(bold("\nNumber of blocks: "), green(length(x)), "\n")

  # Block details
  for (i in seq_along(x)) {
    block_name <- names(x)[i]
    block <- x[[i]]

    # Block header
    cat(bold(blue("\n┌─ Block ")), bold(blue(i)),
        if (!is.null(block_name)) bold(blue(paste0(" (", block_name, ")"))) else "",
        bold(blue(" ─────────────────\n")))

    # Data dimensions
    cat("│ ", bold("Dimensions:"),
        green(paste0(nrow(block$x), " × ", ncol(block$x))), "\n")

    # Design variables
    design_vars <- names(block$design)
    cat("│ ", bold("Design Variables:"),
        green(paste(design_vars, collapse=", ")), "\n")

    # Sample of unique values for each design variable
    cat("│ ", bold("Design Structure:"), "\n")
    for (var in design_vars) {
      unique_vals <- unique(block$design[[var]])
      n_unique <- length(unique_vals)
      sample_vals <- if (n_unique > 5) {
        paste0(paste(head(unique_vals, 3), collapse=", "),
               "...",
               paste(tail(unique_vals, 2), collapse=", "))
      } else {
        paste(unique_vals, collapse=", ")
      }
      cat("│   ", white("•"), " ", var, ": ",
          green(n_unique), " levels (", sample_vals, ")\n", sep="")
    }

    # Column design if present
    if (!is.null(block$column_design)) {
      cat("│ ", bold("Column Design:"), "Present\n")
      col_vars <- names(block$column_design)
      cat("│   Variables: ", green(paste(col_vars, collapse=", ")), "\n")
    }
  }

  # Footer
  cat(bold(blue("\n═══════════════════════\n")))

  invisible(x)
}

#' Print Method for Foldlist Objects
#'
#' Displays a formatted summary of a foldlist object, showing fold structure and held-out information.
#'
#' @param x A foldlist object
#' @param ... Additional arguments passed to print methods
#'
#' @method print foldlist
#' @export
print.foldlist <- function(x, ...) {
  # Get fold information
  foldframe <- attr(x, "foldframe")
  n_folds <- length(x)

  # Header
  cat(crayon::bold(crayon::blue("\n═══ Cross-Validation Folds ═══\n")))

  # Basic information
  cat("\n", crayon::green(paste0("Number of Folds: ", n_folds)), "\n")

  # Print each fold's information
  for (i in seq_len(n_folds)) {
    cat(crayon::bold(paste0("\nFold ", i, ":")), "\n")

    # Get fold details
    fold_info <- x[[i]]

    # Analysis set information
    analysis_size <- nrow(fold_info$analysis[[1]]$design)
    cat("  ", crayon::white("•"), " Analysis Set: ",
        crayon::green(analysis_size), " observations\n", sep="")

    # Assessment set information
    assessment_size <- nrow(fold_info$assessment$design)
    cat("  ", crayon::white("•"), " Assessment Set: ",
        crayon::green(assessment_size), " observations\n", sep="")

    # Held out information if available
    if (!is.null(fold_info$held_out)) {
      cat(crayon::bold("  Held Out Values:"), "\n")
      for (var in names(fold_info$held_out)) {
        val <- fold_info$held_out[[var]]
        cat("    ", crayon::white("◦"), " ", var, ": ",
            crayon::yellow(paste(val, collapse=", ")), "\n", sep="")
      }
    }
  }

  # Footer
  cat(crayon::bold(crayon::blue("\n═════════════════════════════\n")))
  invisible(x)
}
