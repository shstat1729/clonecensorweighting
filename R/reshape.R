#' Create timestamp table
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents a emulated follow up time.
#' @param clone_followup A column name in emulcated clone (i.e. `clones`)
#'   that represents the emulated follow up time in cloned data frame.
#'
#' @returns A data frame with two columns: `tevent` and `ID_t`.
#'   `tevent` represents a timestamp that outcome event can occur based on
#'   observed data. `ID_t` represents an enumerated identifier of each
#'   timestamp, from 1 to n where n represents the number of unique `tevent`
#'   value.
#'
#' @export
#' @examples
create_timestamp_table <- function(clones, clone_followup) {
  timestamps <- sapply(clones, `[[`, i = clone_followup, simplify = FALSE)
  t_events <- sort(unique(unlist(timestamps)))
  res <- dplyr::tibble(tevent = t_events, ID_t = seq_along(t_events))

  res
}


#' Split each observation into multiple subrecords at each time cut
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents a emulated follow up time (corresponding to `clone_followup`
#'   argument) and an event of interest (corresponding to `event` argument).
#' @param clone_followup A column name in emulcated clone (i.e. `clones`)
#'   that represents the emulated follow up time in cloned data frame.
#' @param t_events A vector of timestamp that outcome event can occur based on
#'   observed data.
#' @param event A variable name of an event of interest. The variable should
#'   exists in each data frame that is an element of `clones` argument, and
#'   the variable value should be a binary (0 or 1).
#' @param timestamp_start A new variable name to denote start time.
#' @param id A new variable name for a unique observation identifier, to
#'   represents that multiple rows in output data frame is associated with the
#'   same observation.
#'
#' @returns A list of long-form data frames. Each data frame represents each
#'   clone arm. Each row of the long-form data frame represents a subrecord of
#'   each observation associated with each specific time interval. The first
#'   subrecord starts with time 0, and the rows are expanded up to
#'   `clone_followup`, where cut times are determined by `t_events` argument.
#'
#' @export
#' @examples
split_at_timestamp <- function(
  clones,
  clone_followup,
  t_events,
  event,
  timestamp_start = "Tstart",
  id = "ID"
) {
  arms <- names(clones)

  res <- vector("list", length = length(arms))
  names(res) <- arms
  for (arm in arms) {
    res[[arm]] <-
      clones[[arm]] |>
      survival::survSplit(
        cut = t_events,
        end = clone_followup,
        start = timestamp_start,
        event = event,
        id = id
      )
  }

  res
}


#' Create training data for censoring probability estimation
#'
#' @param clones A list of data frame. Each element of the list represents
#'   each treatment arm. This version of clones must contain a column that
#'   represents an emulated follow up time (corresponding to `clone_followup`
#'   argument), an emulated outcome (correspodning to `clone_outcome`), and a
#'   binary indicator variable that represents whether the observation violates
#'   arm's policy or not (corresponding to `clone_censoring` argument).
#' @param clone_followup A column name that represents the emulated follow up
#'   time in each arm of clones. The variable should exists in each element
#'   data frame of `clones` argument.
#' @param clone_outcome A column name that represents the emulated outcome
#'   in each arm of clones. The variable should exists in each element
#'   data frame of `clones` argument, and the variable value should be binary
#'   (0 or 1).
#' @param clone_censoring A column name that represent whether the observation
#'   violates arm's policy or not. The variable should exists in each element
#'   data frame of `clones` argument, and the variable value should be binary
#'   (0 or 1).
#' @param col_ids A vector of column names that a combination of their values
#'   uniquely identifies each observation.
#' @param timestamp_start A new variable name to denote start time of each
#'   subrecord of observations in a long-form data.
#' @param id A new variable name for a unique observation identifier, to
#'   represents that multiple rows in output data frame is associated with the
#'   same observation.
#' @param timestamp_stop A new variable name to denote end time of each
#'   subrecord of observations in a long-form data.
#'
#' @returns A list of long-form data frames. Each data frame represents each
#'   clone arm. Each row of the long-form data frame represents a subrecord of
#'   each observation associated with each specific time interval. The first
#'   subrecord starts with time 0, and the rows are expanded up to
#'   `clone_followup`, where cut times are determined by `t_events` argument.
#'
#' @export
#' @examples
create_final_data <- function(
  clones,
  clone_followup,
  clone_outcome,
  clone_censoring,
  col_ids,
  timestamp_start = "Tstart",
  id = "ID",
  timestamp_stop = "Tstop"
) {
  df_timestamp <- create_timestamp_table(clones, clone_followup)

  clones_splitted_by_outcome <- split_at_timestamp(
    clones,
    clone_followup,
    df_timestamp$tevent,
    clone_outcome,
    timestamp_start,
    id
  )

  clones_splitted_by_censoring <- split_at_timestamp(
    clones,
    clone_followup,
    df_timestamp$tevent,
    clone_censoring,
    timestamp_start,
    id
  )

  # merge two tables and create column to represent end of timestamp
  df_timestamp_with_time_zero <-
    df_timestamp |>
    dplyr::bind_rows(dplyr::tibble(tevent = 0, ID_t = 0))

  n_clones <- length(clones)
  arms <- names(clones)
  res <- vector("list", length = n_clones)
  names(res) <- arms

  for (i in seq_len(n_clones)) {
    x <-
      clones_splitted_by_outcome[[i]] |>
      dplyr::select(!{{ clone_censoring }})

    y <-
      clones_splitted_by_censoring[[i]] |>
      dplyr::select(
        dplyr::all_of(col_ids),
        {{ clone_followup }},
        {{ clone_censoring }}
      )

    res[[i]] <-
      dplyr::inner_join(
        x,
        y,
        by = dplyr::join_by({{ col_ids }}, {{ clone_followup }})
      ) |>
      dplyr::mutate({{ timestamp_stop }} := .data[[clone_followup]])

    # Merge with timestamp table
    res[[i]] <-
      res[[i]] |>
      dplyr::left_join(
        df_timestamp_with_time_zero,
        by = dplyr::join_by({{ timestamp_start }} == tevent)
      )
  }

  res
}
