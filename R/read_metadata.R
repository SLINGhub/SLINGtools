#' Imports metadata provided by an MSOrganizer EXCEL template
#'
#' @param filename File path of the MSOrganizer EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#'
#' @return A list of tibbles with different metadata
#' @export
#'
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom rlang .data
#' @importFrom tidyselect vars_select_helpers
#' @importFrom dplyr select mutate filter group_by row_number
#'
#' @examples
#' library(SLINGtools)
#'
#' data_file_path <- system.file("extdata",
#'   "Testdata_Lipidomics_MHQuant_Detailed.csv", package = "SLINGtools")
#' d <- read_MassHunterCSV(data_file_path)
#' d
#'
#'
import_MSOrganizerXLM <- function(filename, trim_ws = TRUE){
  d_annot <- list()
  d_annot$annot_analyses <- readxl::read_excel(filename, sheet = "Sample_Annot") |>
    dplyr::mutate(
      VALID_ANALYSIS = TRUE,
      BATCH_ID = as.character(.data$BATCH_ID),
      RUN_ID = dplyr::row_number()) |>
    dplyr::select(
      "RUN_ID",
      ANALYSIS_ID = "Sample_Name",
      DATAFILE_NAME = "Sample_Name",
      QC_TYPE = "Sample_Type",
      SAMPLE_AMOUNT =	"Sample_Amount",
      SAMPLE_AMOUNT_UNIT = "Sample_Amount_Unit",
      ISTD_VOL ="ISTD_Mixture_Volume_[uL]",
      "BATCH_ID",
      "VALID_ANALYSIS"
    ) |>
    dplyr::group_by(.data$BATCH_ID) |>
    dplyr::mutate(BATCH_NO = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  d_annot$annot_features <- readxl::read_excel(filename, sheet = "Transition_Name_Annot", trim_ws = TRUE)|>
    dplyr::mutate(
      FEATURE_ID = stringr::str_squish(.data$Transition_Name),
      FEATURE_NAME = stringr::str_squish(.data$Transition_Name),
      NORM_ISTD_FEATURE_NAME	= stringr::str_squish(.data$Transition_Name_ISTD),
      QUANT_ISTD_FEATURE_NAME = stringr::str_squish(.data$Transition_Name_ISTD),
      isISTD = (.data$FEATURE_ID == .data$Transition_Name_ISTD),
      FEATURE_RESPONSE_FACTOR	= 1,
      isQUANTIFIER = TRUE,
      isINTEGRATED = TRUE,
      REMARKS = NA_character_) |>
    dplyr::select(
      "FEATURE_ID",
      "FEATURE_NAME",
      "isISTD",
      "NORM_ISTD_FEATURE_NAME",
      "QUANT_ISTD_FEATURE_NAME",
      "FEATURE_RESPONSE_FACTOR",
      "isQUANTIFIER",
      "isINTEGRATED",
      "REMARKS")

  #ToDo: Merged cell in template
  annot_istd <- readxl::read_excel(filename,
                           sheet = "ISTD_Annot",
                           skip = 2,
                           trim_ws = TRUE,
                           .name_repair = ~ ifelse(nzchar(.x), .x, LETTERS[seq_along(.x)]))
  names(annot_istd)[1] <- "Transition_Name_ISTD"

  d_annot$annot_istd <- annot_istd |>
    dplyr::mutate(ISTD_COMPOUND_NAME = NA_character_) |>
    dplyr::select(
      QUANT_ISTD_FEATURE_NAME = "Transition_Name_ISTD",
      ISTD_CONC_nM = "ISTD_Conc_[nM]")

  d_annot$annot_responsecurves <- readxl::read_excel(filename, sheet = "Dilution_Annot") |>
    dplyr::select(
      ANALYSIS_ID = "Sample_Name",
      RQC_SERIES_ID = "Dilution_Batch_Name",
      RELATIVE_SAMPLE_AMOUNT = "Relative_Sample_Amount_[%]",
      INJECTION_VOL = "Injection_Volume_[uL]") |>
    dplyr::mutate(
      ANALYSIS_ID = as.character(.data$ANALYSIS_ID),
      RQC_SERIES_ID = as.character(.data$RQC_SERIES_ID),
      RELATIVE_SAMPLE_AMOUNT = .data$RELATIVE_SAMPLE_AMOUNT/100)

  return(d_annot)
}
