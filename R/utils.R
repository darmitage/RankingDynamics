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
#' @description Function to convert an S x T + 1 dataframe of
#' scores into a ranked list. It assumes the data frame's first
#' column is a list of entities and subsequent columns are
#' rankings at each census. Users can either leave N0 empty
#' (which will set N0 as the number of nonzero scores at first census)
#' or specify it (e.g., setting it to 100 will make the list a "top 100" list).
#'
#' @param wide_data A data.frame or a list of data.frames where the first column
#' are species or other identifiers and subsequent columns are their abundances
#' (or more generally scores) at various census points.
#' @param N0 The selected ranking list size. If unspecified, it will calculate
#' N0 as the size of the non-NA entries at the first census (column 2).
#'
#' @return Returns a ranked dataframe or a list of ranked dataframes.
#'
#' @import dplyr
#' @export

get_ranked_data <- function(wide_data, N0 = NULL) {
  # Helper function to process a single dataframe
  process_dataframe <- function(df, N0) {
    # Ensure input is valid
    if (!is.data.frame(df)) stop("Each element in the input list must be a data frame.")
    if (ncol(df) < 2) stop("Each data frame must have at least two columns.")

    # Set default value for N0 if not provided
    if (is.null(N0)) {
      N0 <- sum(!is.na(df[[2]])) # Count non-NA values in the second column
    }

    # Print the value of N0 for debugging
    message("N0 = ", N0)

    # Rank and filter the data
    ranked_data <- df %>%
      mutate(across(-1, as.numeric)) %>%  # Convert all columns (except the first) to numeric
      mutate(across(-1, ~ rank(-., na.last = "keep", ties.method = "random"))) %>% # Rank in descending order
      mutate(across(2, ~ ifelse(. > N0, NA, .))) # Set values above N0 in the second column to NA

    # Remove rows where all ranked columns (except the first) are NA
    filtered_data <- ranked_data[rowSums(!is.na(ranked_data[, -1])) > 0, ]

    return(filtered_data)
  }

  # Check if input is a single dataframe or a list of dataframes
  if (is.data.frame(wide_data)) {
    # Process a single dataframe
    return(process_dataframe(wide_data, N0))
  } else if (is.list(wide_data)) {
    # Process a list of dataframes
    return(lapply(wide_data, process_dataframe, N0 = N0))
  } else {
    stop("Input must be either a data frame or a list of data frames.")
  }
}
