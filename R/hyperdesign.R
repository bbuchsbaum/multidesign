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



#' Function to create a hyperdesign from a data.frame (or tibble)
#'
#' This function takes a data.frame, the names of the design variables, the names of the X variables,
#' and the name of the splitting variable, and returns a hyperdesign object.
#'
#' @importFrom dplyr select everything
#' @importFrom dplyr arrange group_by nest_by
#' @importFrom purrr map
#'
#' @param data A data.frame.
#' @param design_vars a character vector of design variables to extract
#' @param x_vars a character vector of X design variables to extract
#' @param split_var the splitting variable.
#'
#' @return A hyperdesign object.
#'
#' @export
#'
#' @examples
#' # Create a sample tibble
#' sample_tibble <- tibble(
#'   felab = rep(1:2, each = 3),
#'   attention = rep(c("DA", "FA"), times = 3),
#'   basis = rep(c("basis01", "basis02", "basis03"), times = 2),
#'   subject = rep(1001:1002, each = 3),
#'   `1` = rnorm(6),
#'   `2` = rnorm(6),
#'   `3` = rnorm(6)
#' )
#'
#' # Apply the function
#' hd <- df_to_hyperdesign(sample_tibble, c("felab", "attention", "basis"), setdiff(colnames(sample_tibble), c("felab", "attention", "basis", "subject")), "subject")
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



#' Construct a `hyperdesign` object
#'
#' A collection of multivariate datasets (`multidesign` instances) that share a set of common design variables.
#' The class can be used to capture multiblock data, where one wants to model multiple related matrices.
#'
#' @param x a list of `multidesign` instances
#' @param block_names the names of each block
#' @export
#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#'
#' hd <- hyperdesign(list(d1,d2,d3))
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

#' @importFrom multivarious init_transform
#' @return a transformed `hyperdesign` object, with pre-processors set as attribute named `preproc`
#' @export
#' @rdname init_transform
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


#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#'
#' hd <- hyperdesign(list(d1,d2,d3))
#' folds <- fold_over(hd, run)
#' @export
#' @importFrom deflist deflist
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_rows
#' @rdname fold_over
#' folds <- fold_over(hd, run, test_condition = list(task = "a"))
#' @export
#' @importFrom deflist deflist
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_rows
#' @rdname fold_over
fold_over.hyperdesign <- function(x, ..., test_condition = list()) {

  # Validate test_condition if it is non-empty
  if (length(test_condition) > 0) {
    validate_test_condition <- function(x, test_condition) {
      all_design <- do.call(rbind, lapply(x, function(block) block$design))
      for (factor_name in names(test_condition)) {
        if (!factor_name %in% colnames(all_design)) {
          stop(paste("Factor", factor_name, "not found in the design."))
        }
        if (!any(test_condition[[factor_name]] %in% unique(all_design[[factor_name]]))) {
          stop(paste("Level", test_condition[[factor_name]], "not found for factor", factor_name))
        }
      }
    }

    validate_test_condition(x, test_condition)
  }

  splits <- lapply(seq_along(x), function(i) {
    d <- x[[i]]
    split_indices(d, ...) %>% mutate(.block=i)
  })

  foldframe <- splits %>% bind_rows() %>% mutate(.fold=1:n())
  lens <- purrr::map(splits, function(x) nrow(x))
  tlen <- sum(unlist(lens))

  extract <- function(i) {
    block <- foldframe[[".block"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    testdat <- multidesign(x[[block]]$x[ind,], x[[block]]$design[ind,], x[[block]]$column_design)

    traindat <- hyperdesign(lapply(seq_along(x), function(j) {
      if (j == block) {
        if (length(ind) == nrow(x[[j]]$x)) {
          stop(paste("number of test data rows in block ", j, " is equal to number of rows in block.\n                     No training data available for current fold! Check distribution of splitting variables over blocks."))
        }
        multidesign(x[[block]]$x[-ind,], x[[block]]$design[-ind,], x[[block]]$column_design)
      } else {
        x[[j]]
      }
    }))

    # Retrieve the unique values of the held-out split
    held_out_values <- x[[block]]$design[ind, , drop = FALSE] %>%
      distinct() %>%
      slice(1) %>%
      as.list()

    # Retrieve the block indices for the test set
    block_indices <- block_index_mat(x, byrow=TRUE)[block, , drop = FALSE]

    list(analysis = traindat,
         assessment = testdat,
         held_out = held_out_values,
         block_indices = block_indices)
  }

  condition_met <- function(i) {
    block <- foldframe[[".block"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    # Check if the test data meets the specified condition
    test_design <- x[[block]]$design[ind,]
    for (factor_name in names(test_condition)) {
      if (!any(test_design[[factor_name]] %in% test_condition[[factor_name]])) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  indices <- which(sapply(seq_len(tlen), condition_met))
  ret <- deflist(function(i) extract(indices[i]), len=length(indices))

  names(ret) <- paste0("fold_", seq_along(ret))
  class(ret) <- c("foldlist", class(ret))
  attr(ret, "foldframe") <- foldframe
  ret

  ## split over the variable, get global indices
  ## verify that variable exists in each design
  ## verify that all blocks have a minimum observation size
  ## create a `deflist` ,wher each element has an "analysis" and assessment hyperdesign
}

# fold_over.hyperdesign <- function(x, ...) {
#
#   splits <- lapply(seq_along(x), function(i) {
#     d <- x[[i]]
#     split_indices(d, ...) %>% mutate(.block=i)
#   })
#
#
#   foldframe <- splits %>% bind_rows() %>% mutate(.fold=1:n())
#
#   lens <- purrr::map(splits, function(x) nrow(x))
#
#   tlen <- sum(unlist(lens))
#
#   extract <- function(i) {
#     block <- foldframe[[".block"]][i]
#     ind <- unlist(foldframe[["indices"]][[i]])
#
#     testdat <- multidesign(x[[block]]$x[ind,], x[[block]]$design[ind,], x[[block]]$column_design)
#
#     ## all blocks except
#     traindat <- hyperdesign(lapply(seq_along(x), function(j) {
#       if (j == block) {
#         if (length(ind) == nrow(x[[j]]$x)) {
#           stop(paste("number of test data rows in block ", j, " is equal to number of rows in block.
#                      No training data available for current fold! Check distribution of splitting variables over blocks."))
#         }
#         multidesign(x[[block]]$x[-ind,], x[[block]]$design[-ind,], x[[block]]$column_design)
#       } else{
#         x[[j]]
#       }
#     }))
#
#     list(analysis=traindat,
#          assessment=testdat)
#
#   }
#
#
#   ret <- deflist(extract, len=tlen)
#   names(ret) <- paste0("fold_", 1:length(ret))
#   class(ret) <- c("foldlist", class(ret))
#   attr(ret, "foldframe") <- foldframe
#   ret
#
#   ## split over the variable, get global indices
#   ## verify that variable exists in each design
#   ## verify that all blocks have a minimum observation size
#   ## create a `deflist` ,wher each element has an "analysis" and assessment hyperdesign
#
# }



#' @export
xdata.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, xdata)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    xdata(x[[block]])
  }
}

#' @export
design.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, design)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    design(x[[block]])
  }
}

#' @export
column_design.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, design)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, c(1, length(x)))
    column_design(x[[block]])
  }
}

#' @export
subset.hyperdesign <- function(x, fexpr) {
  out <- lapply(hd, function(d) {
    subset(d, !!rlang::enquo(fexpr))
  })

  rem <- unlist(purrr::map(out, is.null))
  if (sum(rem) == length(hd)) {
    stop("subset expression does not match any rows in hyperdesign `x`")
  }

  onam <- names(x)[!rem]
  hyperdesign(out[!rem], onam)
}



#' @export
print.hyperdesign <- function(x) {
  cat("a hyperdesign object. \n")
  cat("with ", length(x), "blocks.", "\n")
  cat("common design variables: ", attr(x, "common_vars"), "\n")
  print(attr(x, "hdes"))
  cat("\n")

}
