#### Utility Functions for preparing and permuting ranked lists

#' @title Permute a ranked list dataframe into N dataframes where ties are randomized
#'
#' @description Since we need to remove ties to calculate ranking metrics, we can attempt to
#' solve the issue of ties (common in ecological rankings), by permuting scores
#' a bunch of times, randomly assigning arbitrarily small values to elements that are tied. This allows us to assess
#' the sensitivity of the results on our ordering of entries.
#'
#' @param wide_data A data.frame object where the first column are species or other
#' identifiers and subsequent columns are their abundances (or more generally scores)
#' at various census points. Thus, the dataframe should be S (species no) x T (census times) + 1
#' @param N The selected ranking list size. If unspecified, it will just calculate N0 as the size of the
#' non-NA entries at the first census (column 2).
#' @param epsilon The small values added randomly to tied scores in order to remove ties. If not provided,
#' it is calculated as one-hundredth of the smallest non-NA value in wide_data[,-1].
#' @return Returns a list of length N of ranked_lists, each one a new permutation where
#' ties are removed and randomly reordered.
#'
#' @import dplyr
#' @export
#'
permute_dataframe <- function(wide_data, N, epsilon = NULL) {
  # Validate inputs
  if (!is.data.frame(wide_data)) stop("Input 'wide_data' must be a data frame.")
  if (!is.numeric(N) || N <= 0 || N != as.integer(N)) stop("'N' must be a positive integer.")

  # Calculate epsilon if not provided
  if (is.null(epsilon)) {
    non_na_values <- wide_data[,-1][!is.na(wide_data[,-1])]
    if (length(non_na_values) == 0) stop("No non-NA values in the dataset to calculate epsilon.")
    epsilon <- min(non_na_values) / 100
    message("Epsilon automatically set to: ", epsilon)
  } else if (!is.numeric(epsilon) || epsilon <= 0) {
    stop("'epsilon' must be a positive number.")
  }

  # Create a list to store the permuted dataframes
  permuted_list <- vector("list", N)

  # Perform permutations
  for (i in 1:N) {
    # Make a copy of the dataframe for this permutation
    df_permuted <- wide_data

    # Loop through each column, except the 'entity' column
    for (col in names(wide_data)[-1]) {
      # Identify duplicated values in the column (ignoring NAs)
      dup_indices <- which(duplicated(df_permuted[[col]]) & !is.na(df_permuted[[col]]))

      if (length(dup_indices) > 0) {
        # Generate small random values to add to the duplicated rows
        random_values <- runif(length(dup_indices), min = 0, max = epsilon)

        # Add the random values to break ties
        df_permuted[[col]][dup_indices] <- df_permuted[[col]][dup_indices] + random_values
      }
    }

    # Store the permuted dataframe in the list
    permuted_list[[i]] <- df_permuted
  }

  return(permuted_list)
}

#' @title Create a ranked list from a table of scores
#'
#' @description Convert an S x T+1 data.frame of scores (species x time points)
#' into a ranked list of size N0. If N0 is not provided, it's set to the minimum number
#' of non-NA entries across all time points, ensuring consistency.
#'
#' @param wide_data A data.frame or list of data.frames with species in column 1,
#' and abundance/score values in columns 2 onward.
#' @param N0 Integer, number of ranks to retain. If NULL, uses min non-NA count per column.
#'
#' @return A ranked data.frame or list of ranked data.frames with only top N0 per column retained.
#' @import dplyr
#' @export

get_ranked_data <- function(wide_data, N0 = NULL) {

  process_dataframe <- function(df, N0) {
    if (!is.data.frame(df)) stop("Each element must be a data.frame.")
    if (ncol(df) < 2) stop("Each data.frame must have at least two columns.")

    # Convert scores to numeric
    df <- df %>%
      mutate(across(-1, as.numeric))

    # Determine N0 if not provided: minimum number of non-NA entries across time points
    if (is.null(N0)) {
      N0 <- min(colSums(!is.na(df[, -1])))
    }

    message("Using N0 = ", N0)

    # Rank each time point and mask values beyond N0
    ranked <- df %>%
      mutate(across(-1, ~ rank(-., na.last = "keep", ties.method = "random"))) %>%
      mutate(across(-1, ~ ifelse(. > N0, NA, .)))

    # Remove rows that are NA across all time points
    ranked <- ranked[rowSums(!is.na(ranked[, -1])) > 0, ]

    return(ranked)
  }

  # Dispatch based on input type
  if (is.data.frame(wide_data)) {
    return(process_dataframe(wide_data, N0))
  } else if (is.list(wide_data)) {
    return(lapply(wide_data, process_dataframe, N0 = N0))
  } else {
    stop("Input must be a data.frame or list of data.frames.")
  }
}
