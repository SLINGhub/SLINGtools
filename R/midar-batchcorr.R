#' Combines a list of MidarExperiments into one
#'
#' @param data MidarExperiment object
#' @param qc_types QC types used for batch correction
#' @param center_by_fun Function used to center. Default is "median".
#' @export
#'
#' @importFrom glue glue
#' @importFrom openxlsx write.xlsx
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion
#'

batch_corr_center <- function(data, qc_types, center_by_fun = "median"){

  ds <- data@dataset
  ds <- ds %>%
    group_by(.data$FEATURE_NAME,  .data$BATCH_ID) %>%
    mutate(CONC_FINAL = .data$Concentration/median(.data$Concentration[.data$QC_TYPE %in% qc_types], na.rm = TRUE)) %>%
    ungroup()


  ds <- ds %>%
    group_by(.data$FEATURE_NAME) %>%
    mutate(CONC_FINAL =  CONC_FINAL * median(Concentration[.data$QC_TYPE %in% qc_types], na.rm = TRUE))

  data@dataset <- ds
  data
}


#' Performs drift correction based on LOESS smoothing
#'
#' @param data MidarExperiment object
#' @param qc_types QC types used for drift correction
#' @param smooth_by_batch Correct each batch separately (Default is TRUE)
#' @param log2_transform log2 transform data during correction (Default is TRUE)
#' @param loess_span Loess span width (default is 0.75)
#' @param limit_to_species Apply correction only to species matching (RegEx)
#' @param apply_only_if_smpl_cv_lower Apply correction only if sample CV decreases (default is FALSE)
#' @param Minimum CV difference when using apply_only_if_smpl_cv_lower
#'
#' @export
#'


drift_corr_loess <- function(data, qc_types, smooth_by_batch = TRUE, log2_transform = TRUE, loess_span = 0.75,
                             limit_to_species = NULL, apply_only_if_smpl_cv_lower = FALSE, treshold_smpl_cv_delta = 0){
 # browser()
  get_loess <- function(d, qc_types,loess_span) {

    tryCatch({
      res <- stats::loess(y ~ x,
                          span = loess_span, family = "gaussian", degree = 2, normalize=FALSE, iterations=4,
                          data = d[d$QC_TYPE %in% qc_types, ]) %>%
        stats::predict(tibble::tibble(x = seq(min(d$x), max(d$x), 1))) %>% as.numeric()
      res
    },
    error = function(e) {
      return(rep(NA_real_, length(d$x)))
    })
  }


  if(is.null(limit_to_species))
    ds <- data@dataset
  else
    ds <- data@dataset %>% dplyr::filter(stringr::str_detect(FEATURE_NAME, limit_to_species))


  ds$x <- ds$RUN_ID
  ds$y <- ds$Concentration
  if(log2_transform) ds$y <- log2(ds$y)

  d <- ds %>%
    group_by(.data$FEATURE_NAME, .data$BATCH_ID) %>%
    nest() %>%
    mutate(Y_PREDICTED = purrr::map(data, \(x) get_loess(x, qc_types, loess_span))) %>%
    unnest(cols = c(data, Y_PREDICTED))


  if(log2_transform){
    d <- d %>%
      group_by(FEATURE_NAME, .data$BATCH_ID) %>%
      mutate(Y_PREDICTED = Y_PREDICTED - median(Y_PREDICTED, na.rm = TRUE),
             Y_ADJ = 2^(.data$y - Y_PREDICTED))
  } else {
    d <- d %>%
      group_by(FEATURE_NAME, .data$BATCH_ID) %>%
      mutate(Y_PREDICTED = Y_PREDICTED / median(Y_PREDICTED, na.rm = TRUE),
             Y_ADJ = .data$y / Y_PREDICTED)
  }

  if(!apply_only_if_smpl_cv_lower)  treshold_smpl_cv_delta = Inf

  d <- d %>%
    group_by(FEATURE_NAME, .data$BATCH_ID) %>%
    mutate(CV_RAW_SPL = sd(Concentration[QC_TYPE == "SPL"], na.rm = TRUE)/mean(Concentration[QC_TYPE == "SPL"], na.rm = TRUE) *100,
           CV_ADJ_SPL = sd(Y_ADJ[QC_TYPE == "SPL"], na.rm = TRUE)/mean(Y_ADJ[QC_TYPE == "SPL"], na.rm = TRUE) *100) %>%
    mutate(
      DRIFT_CORRECTED = (CV_RAW_SPL - CV_ADJ_SPL) >  treshold_smpl_cv_delta,
      Y_FINAL = if_else(apply_only_if_smpl_cv_lower & DRIFT_CORRECTED, Y_ADJ, Concentration))


  data@dataset <- data@dataset %>% dplyr::left_join(d %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", CURVE_PREDICTED = "Y_PREDICTED", CONC_DRIFT_ADJ = "Y_ADJ", "CV_RAW_SPL", "CV_ADJ_SPL", "DRIFT_CORRECTED", CONC_FINAL = "Y_FINAL"))
  data
}
