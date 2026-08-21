reference_masked_summary <- function(x, cells, groups, fun) {
  levels <- unique(groups)
  out_x <- matrix(NA_real_, nrow = length(levels), ncol = ncol(x))
  out_cells <- matrix(FALSE, nrow = length(levels), ncol = ncol(x))
  for (i in seq_along(levels)) {
    rows <- which(groups == levels[[i]])
    for (j in seq_len(ncol(x))) {
      observed <- cells[rows, j]
      if (any(observed)) {
        out_x[i, j] <- fun(x[rows[observed], j])
        out_cells[i, j] <- TRUE
      }
    }
  }
  list(x = out_x, cells = out_cells)
}

test_that("masked summaries require an explicit aggregation rule", {
  x <- matrix(seq_len(24), nrow = 6, ncol = 4)
  groups <- rep(c("A", "B"), each = 3)
  cells <- matrix(seq_len(length(x)) %% 3L != 0L, nrow = nrow(x))
  md <- multidesign(x, data.frame(group = groups), cells = cells)

  expect_error(summarize_by(md, group), "requires an explicit `aggregate`")

  result <- summarize_by(md, group, aggregate = "mean")
  reference <- reference_masked_summary(x, cells, groups, mean)
  expect_equal(result$x, reference$x, ignore_attr = TRUE)
  expect_identical(unname(cell_mask(result)), unname(reference$cells))
  expect_equal(design(result)$group, c("A", "B"))
})

test_that("masked summaries represent fully unobserved coordinates explicitly", {
  x <- matrix(seq_len(12), nrow = 4, ncol = 3)
  groups <- c("A", "A", "B", "B")
  cells <- matrix(TRUE, nrow = 4, ncol = 3)
  cells[1:2, 2] <- FALSE
  x[3, 1] <- NA_real_
  md <- multidesign(x, data.frame(group = groups), cells = cells)

  result <- summarize_by(md, group, aggregate = "mean")
  expect_true(is.na(result$x[1, 2]))
  expect_false(cell_mask(result)[1, 2])
  expect_true(is.na(result$x[2, 1]))
  expect_true(cell_mask(result)[2, 1])
})

test_that("custom masked aggregation is checked and unmasked legacy is unchanged", {
  x <- matrix(seq_len(16), nrow = 4)
  groups <- c("A", "A", "B", "B")
  cells <- matrix(TRUE, nrow = 4, ncol = 4)
  cells[1, 1] <- FALSE
  masked <- multidesign(x, data.frame(group = groups), cells = cells)

  maximum <- summarize_by(masked, group, aggregate = max)
  reference <- reference_masked_summary(x, cells, groups, max)
  expect_equal(maximum$x, reference$x, ignore_attr = TRUE)
  expect_identical(unname(cell_mask(maximum)), unname(reference$cells))

  expect_error(
    summarize_by(masked, group, aggregate = function(values) character()),
    "group 1.*coordinate 1.*one numeric"
  )

  unmasked <- multidesign(x, data.frame(group = groups))
  legacy <- summarize_by(unmasked, group)
  expect_equal(legacy$x, rbind(colMeans(x[1:2, ]), colMeans(x[3:4, ])))
  expect_null(cell_mask(legacy))

  explicit <- summarize_by(unmasked, group, aggregate = "mean")
  expect_equal(explicit$x, legacy$x, ignore_attr = TRUE)
  expect_null(cell_mask(explicit))
})

test_that("masked aggregation matches a randomized reference loop", {
  set.seed(20260822)
  for (iteration in seq_len(25L)) {
    n <- sample(4:10, 1L)
    p <- sample(2:6, 1L)
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    cells <- matrix(sample(c(TRUE, FALSE), n * p, replace = TRUE), nrow = n)
    groups <- sample(LETTERS[1:3], n, replace = TRUE)
    groups[seq_len(min(3L, n))] <- LETTERS[seq_len(min(3L, n))]
    md <- multidesign(x, data.frame(group = groups), cells = cells)

    actual <- summarize_by(md, group, aggregate = "mean")
    reference <- reference_masked_summary(x, cells, groups, mean)
    expect_equal(actual$x, reference$x, tolerance = 1e-12, ignore_attr = TRUE)
    expect_identical(unname(cell_mask(actual)), unname(reference$cells))
  }
})
