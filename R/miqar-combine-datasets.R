
#' Combines a list of MidarExperiments into one
#'
#' @param ... MidarExperiment objects
#' @param run_order_as_list File name and path of the Excel file
#' @export
#'
#' @importFrom glue glue
#' @importFrom openxlsx write.xlsx
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion
#'
#'


combine_experiments <- function(..., ordered_by_runsequence){

  exp_list <- list(...)

  mexp <- MidarExperiment()


  mexp@dataset <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@dataset)  |> distinct()
  mexp@annot_analyses <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_analyses) |> mutate(RUN_ID_ANNOT = row_number())
  mexp@annot_istd <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_istd)  |> distinct()
  mexp@annot_features <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_features) |> distinct()

  mexp@dataset <- mexp@dataset %>%
    rename(BATCH_RUN_ID = RUN_ID) %>%
    group_by(FEATURE_NAME) %>%
    mutate(RUN_ID = row_number(), .before = BATCH_RUN_ID) %>%
    ungroup()

  mexp@annot_batch_info <- mexp@annot_analyses %>%
    dplyr::group_by(.data$BATCH_ID) %>%
    dplyr::mutate(BATCH_NO = dplyr::cur_group_id()) %>%
    dplyr::summarise(
      BATCH_ID = .data$BATCH_ID[1],
      BATCH_NO = .data$BATCH_NO[1],
      id_batch_start = dplyr::first(.data$RUN_ID_ANNOT),
      id_batch_end = dplyr::last(.data$RUN_ID_ANNOT)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)

  mexp

}









