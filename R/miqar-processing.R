#' get_conc_unit
#'
#' @param sample_amount_unit MidarExperiment object
#' @return string with concentration unit
#' @noRd

get_conc_unit <- function(sample_amount_unit){
  if (length(unique(sample_amount_unit)) > 1)
    conc_unit <- "pmol/sample amount unit (multiple units)"
  else if (sample_amount_unit[1] == "uL" | sample_amount_unit[1] == "\U003BCL")
    conc_unit <- "\U003BCmol/L"
  else
    conc_unit <- glue::glue("pmol/{sample_amount_unit}")
  conc_unit
}


#' normalizeByISTD
#'
#' @param data MidarExperiment object
setGeneric("normalizeByISTD", function(data) standardGeneric("normalizeByISTD"))



#' Normalize Intensities with corresponding ISTD Intensities
#'
#' @param data MidarExperiment object
#' @return MidarExperiment object
#'
#' @export
#'
#' @importFrom glue glue

setMethod("normalizeByISTD", signature = "MidarExperiment", function(data) {
  if(nrow(data@annot_features) < 1) stop("ISTD map is missing...please import transition annotations.")

  d_temp <- data@dataset %>%
    dplyr::full_join(data@annot_features, by = c("FEATURE_NAME"))

  d_temp <- d_temp  %>%
    dplyr::group_by(.data$NORM_ISTD_FEATURE_NAME, .data$ANALYSIS_ID) %>%
    dplyr::mutate(normIntensity = .data$Intensity/.data$Intensity[.data$isISTD]) %>%
    dplyr::ungroup()

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", "isISTD", "normIntensity"), by = c("ANALYSIS_ID", "FEATURE_NAME"))

  n_features <- length(unique(data@annot_features$FEATURE_NAME))
  n_ISTDs <- length(unique(data@annot_features$NORM_ISTD_FEATURE_NAME))
  print(glue::glue("{n_features} features normalized with {n_ISTDs} ISTDs. \n"))

  data
})

#' quantitateByISTD
#'
#' @param data MidarExperiment object
setGeneric("quantitateByISTD", function(data) standardGeneric("quantitateByISTD"))


#' Quantitate using sample and spiked ISTD amounts
#'
#' @param data MidarExperiment object
#' @return MidarExperiment object
#'
#' @export
#'
#' @importFrom glue glue
#' @importFrom rlang .data

setMethod("quantitateByISTD", signature = "MidarExperiment", function(data) {
  if(nrow(data@annot_istd) < 1) stop("ISTD concetrations are missing...please import them first.")

  if(!(c("normIntensity") %in% names(data@dataset))) stop("Data needs first to be ISTD normalized. Please apply function 'normalizeByISTD' first.")

  d_temp <- data@dataset %>%
    dplyr::full_join(data@annot_features, by = c("FEATURE_NAME")) %>%
    dplyr::full_join(data@annot_istd, by = c("QUANT_ISTD_FEATURE_NAME"))

  d_temp <- d_temp %>% mutate(pmol_total = (.data$normIntensity)*(.data$ISTD_VOL*(.data$ISTD_CONC_nM))/1000)
  d_temp <- d_temp %>% mutate(Concentration = .data$pmol_total/.data$SAMPLE_AMOUNT)

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", "pmol_total", "Concentration"), by = c("ANALYSIS_ID", "FEATURE_NAME"))

  n_features <- length(unique(data@annot_features$FEATURE_NAME))
  n_ISTDs <- length(unique(data@annot_features$NORM_ISTD_FEATURE_NAME))

  conc_unit <- get_conc_unit(data@annot_analyses$SAMPLE_AMOUNT_UNIT)

  print(glue::glue("{n_features} compounds quantitated in {nrow(data@annot_analyses)} samples using {n_ISTDs} spiked ISTDs.
                   Concentration unit: [{conc_unit}]"))

  data
})



#' exportWideCSV generic
#'
#' @param data MidarExperiment object
#' @param variable Variable to be exported
#' @param filename File name with path of exported CSV file

setGeneric("exportWideCSV", function(data, variable, filename) standardGeneric("exportWideCSV"))

#' Export any parameter to a wide-format table
#'
#' @param data MidarExperiment object
#' @param variable Variable to be exported
#' @param filename File name with path of exported CSV file
#'
#' @export
#'
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#'
#'
setMethod("exportWideCSV", signature = "MidarExperiment", function(data, variable, filename) {

  var <- dplyr::sym(variable)

  if (!(variable %in% names(data@dataset))) stop("Variable '", variable,  "' does not (yet) exist in dataset")

  ds <- data@dataset |>
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "AcqTimeStamp", "FEATURE_NAME", !!var) %>%
    tidyr::pivot_wider(names_from = .data$FEATURE_NAME, values_from = !!var)

  readr::write_csv(ds, file = filename, num_threads = 4, col_names = TRUE)
  invisible(ds)

})

#' calcQC generic
#'
#' @param data MidarExperiment objec

setGeneric("calcQC", function(data) standardGeneric("calcQC"))

#' Calculate QC parameters for each feature
#'
#' @param data MidarExperiment object
#'
#' @export
#'
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider nest unnest
#' @importFrom rlang .data
#' @importFrom purrr map
#' @importFrom broom glance
#' @importFrom dplyr summarise
#'
#'
setMethod("calcQC", signature = "MidarExperiment", function(data) {

  if(!(c("normIntensity") %in% names(data@dataset))) stop("Data needs first to be ISTD normalized. Please apply function 'normalizeByISTD' first.")

  ds1 <- data@dataset %>%
    dplyr::filter(.data$QC_TYPE %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK")) %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::summarise(
      Int_med_PBLK = median(.data$Intensity[.data$QC_TYPE == "PBLK"], na.rm = TRUE),
      Int_med_SPL = median(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE),
      Int_med_BQC = median(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE),
      Int_med_TQC = median(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE),
      Int_med_NIST = median(.data$Intensity[.data$QC_TYPE == "NIST"], na.rm = TRUE),
      Int_med_LTR = median(.data$Intensity[.data$QC_TYPE == "LTR"], na.rm = TRUE),

      SB_Ratio_Q10 = quantile(.data$Intensity[.data$QC_TYPE == "SPL"], probs  = 0.1, na.rm = TRUE, names = FALSE)/median(.data$Intensity[.data$QC_TYPE == "PBLK"], na.rm = TRUE, names = FALSE),

      Int_CV_TQC = sd(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE) * 100,
      Int_CV_BQC = sd(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE) * 100,
      Int_CV_SPL = sd(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE) * 100,

      normInt_CV_TQC = sd(.data$normIntensity[.data$QC_TYPE == "TQC"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "TQC"], na.rm = TRUE) * 100,
      normInt_CV_BQC = sd(.data$normIntensity[.data$QC_TYPE == "BQC"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "BQC"], na.rm = TRUE) * 100,
      normInt_CV_SPL = sd(.data$normIntensity[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "SPL"], na.rm = TRUE) * 100,
    )

  model <- as.formula("Intensity ~ RELATIVE_SAMPLE_AMOUNT")

  ds2 <- data@dataset %>%
    dplyr::filter(.data$QC_TYPE %in% c("RQC")) %>%
    dplyr::full_join(data@annot_responsecurves, by = "ANALYSIS_ID") %>%
    dplyr::group_by(.data$FEATURE_NAME, .data$RQC_SERIES_ID) %>%
    dplyr::filter(!all(is.na(.data$Intensity))) %>%
    tidyr::nest() %>%
      mutate(
        models = purrr::map(data, function(x) lm(model, data = x, na.action = na.exclude)),
        #mandel = map(data, \(x) DCVtestkit::calculate_mandel(x, "RELATIVE_SAMPLE_AMOUNT", "Intensity")),
        #ppa = map(data, \(x) DCVtestkit::calculate_pra_linear(x, "RELATIVE_SAMPLE_AMOUNT", "Intensity")),
        tidy = purrr::map(.data$models, function(x) broom::glance(x))) %>%
    tidyr::unnest(c("tidy")) %>%
    dplyr::select("FEATURE_NAME", "RQC_SERIES_ID", R2 = "r.squared", Y0 = "sigma") %>%
    tidyr::pivot_wider(names_from = "RQC_SERIES_ID", values_from = c("R2", "Y0"), names_prefix = "RQC_")


  data@d_QC <- ds1 %>% dplyr::left_join(ds2, by = "FEATURE_NAME")
  data

})


#' saveQCinfo generic
#'
#' @param data MidarExperiment object
#' @param filename File name with path of exported CSV file

setGeneric("saveQCinfo", function(data, filename) standardGeneric("saveQCinfo"))

#' Save the QC table to a CSV file
#'
#' @param data MidarExperiment object
#' @param filename File name with path of exported CSV file
#'
#' @export
#'
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#'
#'
setMethod("saveQCinfo", signature = "MidarExperiment", function(data, filename) {

  if (nrow(data@d_QC)== 0) stop("QC info has not yet been calculated. Please apply 'calcQC' first.")

  readr::write_csv(data@d_QC, file = filename, num_threads = 4, col_names = TRUE)
  invisible(data@d_QC)

})

#' getDatasetFilteredQC generic
#'
#' @param data MidarExperiment object
#' @param Intensity_BQC_min Intensity_BQC_min
#' @param CV_BQC_max = CV_BQC_max
#' @param Intensity_TQC_min = Intensity_TQC_min
#' @param CV_TQC_max = CV_TQC_max
#' @param SB_RATIO_min = SB_RATIO_min
#' @param R2_min = R2_min
#' @param RQC_CURVE = RQC_CURVE
#'
setGeneric("getDatasetFilteredQC", function(data, Intensity_BQC_min = NA, CV_BQC_max = NA, Intensity_TQC_min = NA, CV_TQC_max = NA, SB_RATIO_min = NA, R2_min = NA, RQC_CURVE = NA) standardGeneric("getDatasetFilteredQC"))

#' Save the QC table to a CSV file
#'
#' @param data MidarExperiment object
#' @param Intensity_BQC_min Intensity_BQC_min
#' @param CV_BQC_max = CV_BQC_max
#' @param Intensity_TQC_min = Intensity_TQC_min
#' @param CV_TQC_max = CV_TQC_max
#' @param SB_RATIO_min = SB_RATIO_min
#' @param R2_min = R2_min
#' @param RQC_CURVE = RQC_CURVE
#'
#' @export
#'
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @importFrom stats as.formula lm median na.exclude quantile sd
#'
#'
setMethod("getDatasetFilteredQC",
          signature = "MidarExperiment",
          function(data,
                   Intensity_BQC_min = NA,
                   CV_BQC_max = NA,
                   Intensity_TQC_min = NA,
                   CV_TQC_max = NA,
                   SB_RATIO_min = NA,
                   R2_min = NA,
                   RQC_CURVE = 1
                   ){

  if (nrow(data@d_QC)== 0) stop("QC info has not yet been calculated. Please apply 'calcQC' first.")

  if(is.na(Intensity_BQC_min)) Intensity_BQC_min <- -Inf
  if(is.na(CV_BQC_max)) CV_BQC_max <- Inf
  if(is.na(Intensity_TQC_min)) Intensity_TQC_min <- -Inf
  if(is.na(CV_TQC_max)) CV_TQC_max <- Inf
  if(is.na(SB_RATIO_min)) SB_RATIO_min <- 0
  if(is.na(R2_min)) R2_min <- 0

  d_filt <-  data@d_QC %>% filter(.data$Int_med_BQC > Intensity_BQC_min,
                                  .data$Int_med_TQC > Intensity_TQC_min,
                                  .data$normInt_CV_BQC < CV_BQC_max,
                                  .data$normInt_CV_TQC < CV_TQC_max,
                                  .data$SB_Ratio_Q10 > SB_RATIO_min)

  print(glue::glue("{nrow(d_filt)} of {nrow(data@d_QC)} features passed QC filtering."))
  data@dataset_QC_filtered <- data@dataset %>% dplyr::right_join(d_filt|> dplyr::select("FEATURE_NAME"), by = "FEATURE_NAME")
  data
})



