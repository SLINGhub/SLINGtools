#' Read and convert an Agilent MassHunter Quant CSV result file
#'
#' @param file File path of MassHunter Quant CSV file
#' @param silent Suppress messages
#'
#' @return A tibble in the long format
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
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
read_MassHunterCSV <- function(file, silent = FALSE) {

  if(!silent) cat(paste0("Reading '", basename(file), "' ... "), fill = TRUE, sep = " ")
  # if (shiny::isRunning())
  #   incProgress(1 / length(n_datafiles), detail = paste0("", basename(file)))
  #
  # Read Agilent MassHunter Quant Export file (CSV)
  datWide <-
    readr::read_csv(
      file = file,
      col_names = FALSE,
      na = c("#N/A", "NULL"),
      trim_ws = TRUE,
      col_types = readr::cols(.default = "c"),
      locale = readr::locale(encoding = 'ISO-8859-1'), num_threads = 4,progress = TRUE
    )

  # Remove text that is not required and remove dot chars that interfere later with the conversion wide to long
  # ToDo: Convert to tidyverse functions
  datWide[1, ] <-
    lapply(datWide[1, ], function(y)
      gsub(" Results", "", y))
  datWide[1, ] <-
    lapply(datWide[1, ], function(y)
      gsub(" Method", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("\\. ", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("\\.", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("/", "", y))


  datWide <- datWide |> dplyr::add_row(.after = 1)
  datWide[1, ] <- tibble::tibble(A = datWide[1,] |> unlist() |> dplyr::na_if("")) |>  tidyr::fill(A) |> unlist() |> as.list()

  datWide[1, ] <- replace(datWide[1,], stringr::str_detect(string = datWide[1,] , pattern = "AA"),"")
  datWide[2, ] <- replace(datWide[1,], !stringr::str_detect(string = datWide[1,] , pattern = "AA"),"")
  datWide[1, ] <- tibble::tibble(A = datWide[1,] |> unlist() |> dplyr::na_if("")) |>  tidyr::fill(A) |> unlist() |> as.list()

   # Concatenate rows
  datWide[1, ] <- paste(datWide[1, ], datWide[2, ], sep = " ") |> stringr::str_squish() |> as.list()
  datWide <- datWide[-2, ]


  # Fill in transition name in empty columns (different parameters of the same transition)
  # Dealing with exported "Qualifier" results: Adds the corresponding quantifier name (the first column of the group) with the Tag "QUAL" and the transition in front e.g. "Sph d18:0 [QUAL 302.3>266.2]"


  # Remove columns with no column names, as they are undefined (seems to be only Outlier Summary and Quantitation Message Summary)
  datWide <- datWide |> dplyr::select(-where(~ .x[2] == ""))

  # Concatenate rows containing parameters + transitions to the form parameter.sample and parameter.transition headers. This will later be converted with reshape()
  colnames(datWide) <- paste(datWide[2, ], datWide[1, ], sep = "\t")
  datWide <- datWide[-1:-2, ]


  # Rename some known column header names and remove columns that are not needed or not informative.

  new_colnames = c(
    "DataFileName" = "Data File\tSample",
    "SampleName" = "Name\tSample",
    "AcqTimeStamp" = "AcqDate-Time\tSample",
    "SampleType" = "Type\tSample",
    "VialPosition" = "Pos\tSample",
    "InjVolume" = "Vol\tSample",
    "Dilution" = "Dil\tSample",
    "SampleGroup" = "Sample Group\tSample",
    "InstrumentName" = "Instrument Name\tSample",
    "InstrumentType" = "Instrument Type\tSample",
    "AcqMethodFile" = "Acq. Method File\tSample",
    "AcqMethodPath" = "Acq. Method Path\tSample",
    "DataPath" = "Data Path\tSample",
    "QuantitationMessage" = "Quantitation Message\tSample",
    "NameCompound" = "Name\tCompound"
  )

  datWide <- datWide %>% dplyr::rename(dplyr::any_of(new_colnames))

  if("QuantitationMessage" %in% names(datWide)) stop("Field 'Quantitation Message' currently not supported: Please re-export your data in MH without this field.")
  if("NameCompound" %in% names(datWide)) stop("Compound table format is currently not supported. Please re-export your data in MH with compounds as columns.")

  if(! "DataFileName" %in% names(datWide)) stop("Unknown format, or data file is corrupt. Please try re-export from MH.")
  # Remove ".Sample" from remaining sample description headers and remove known unused columns
  datWide <-
    datWide[, !(names(datWide) %in% c("NA\tSample", "Level\tSample", "Sample"))]
  names(datWide) <- sub("\\\tSample", "", names(datWide))

  keep.cols <- names(datWide) %in% c("")
  datWide <- datWide [!keep.cols]



  # Transform wide to long (tidy) table format
  # ------------------------------------------



  datWide <- datWide %>%
    dplyr::mutate(
      dplyr::across(.cols = dplyr::any_of(c("AcqTimeStamp")),
                    .fns = lubridate::mdy_hm),
      dplyr::across(.cols = dplyr::any_of(c("SampleName")),
                    .fns = stringr::str_squish)
    )

  # obtain list with column names of the columns that define transition values (e.g. "RT Cer d16:0/18:0"). Delimuter is currently tab (\t)
  param_transition_names <-
    colnames(datWide[, -1:-tail(grep("\\\t", colnames(datWide), invert =  TRUE), 1)])


  # Obtain long table of all param-transition combinations, split param and compund name and then spread values of different param as columns
  datLong <- datWide %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(param_transition_names), names_pattern = "(.*)\t(.*)$", names_to = c("Param", "Feature")) %>%
    tidyr::pivot_wider(names_from = "Param" ,values_from = "value")

  # Convert types of knows parameters and fields in the data set
  # ------------------------------------------------------------

  datLong <- datLong %>%
    dplyr::mutate(
      dplyr::across(.cols = dplyr::any_of(c("RT", "Area", "Height","FWHM","Width","SN","IntStart","IntEnd",
                                    "Symmetry","InjVolume", "Precursor Ion", "Product Ion", "Collision Energy")),
                    .fns = as.numeric),
      dplyr::across(.cols = dplyr::any_of(c("MI")),
                    .fns = as.logical),
      dplyr::across(.cols = dplyr::any_of(c("Ion Polarity")),
                    .fns = as.factor)
      )

  new_colnames <- c(DataName = "SampleName", PrecursorMZ = "Precursor Ion", ProductMZ = "Product Ion",
                    CollisionEnergy = "Collision Energy", IonPolarity = "Ion Polarity")

  datLong <- datLong %>% dplyr::rename(any_of(new_colnames))

  if(!silent) {
    cat("Imported ", length(unique(datLong$DataFileName)), "samples with ", fill = FALSE)
    cat(length(unique(datLong$Feature)), "transitions ", fill = FALSE)
  }

  datLong
}


#' Reads a defined peak parameter from a Agilent MassHunter Quant CSV result file into a flat wide table
#'
#' @param file File name and path of the MassHunter Quant CSV file
#' @param field Peak parameter (e.g. Area, RT)
#' @param silent Suppress messages
#'
#' @return A tibble in the long format
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers

#'
#' @examples
#' library(SLINGtools)
#'
#' data_file_path <- system.file("extdata",
#'   "Testdata_Lipidomics_MHQuant_Detailed.csv", package = "SLINGtools")
#' d_area <- read_MassHunterCSV_wide(data_file_path, field = "Area")
#' d_area
#'
#'
read_MassHunterCSV_wide <- function(file, field, silent = FALSE) {

  sample_def_cols = c(
    "DataFileName",
    "SampleName",
    "AcqTimeStamp",
    "SampleType",
    "VialPosition",
    "InjVolume",
    "SampleType"
  )
  d <- read_MassHunterCSV(file, silent) %>%
    dplyr::select(tidyselect::any_of(sample_def_cols), .data$Feature, {{field}})

  d %>% tidyr::pivot_wider(names_from = "Feature", values_from = {{field}})
}
