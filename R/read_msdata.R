#' Read and convert an Agilent MassHunter Quant CSV result file
#'
#' @param rawFileName File path of MassHunter Quant CSV file
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
read_MassHunterCSV <- function(rawFileName, silent = FALSE) {

  if(!silent) cat(paste0("Reading '", basename(rawFileName), "' ... "), fill = TRUE, sep = " ")
  # if (shiny::isRunning())
  #   incProgress(1 / length(n_datafiles), detail = paste0("", basename(rawFileName)))
  #
  # Read Agilent MassHunter Quant Export file (CSV)
  datWide <-
    readr::read_csv(
      file = rawFileName,
      col_names = FALSE,
      na = c("#N/A", "NULL"),
      trim_ws = TRUE,
      col_types = readr::cols(.default = "c"),
      locale = readr::locale(encoding = 'ISO-8859-1')
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


  # Fill in transition name in empty columns (different parameters of the same transition)
  # Dealing with exported "Qualifier" results: Adds the corresponding quantifier name (the first column of the group) with the Tag "QUAL" and the transition in front e.g. "Sph d18:0 [QUAL 302.3>266.2]"
  col_name = "Sample"
  isFirstQualifier = TRUE
  quantifier_name =
  for (c in 1:ncol(datWide)) {
    val = datWide[1, c]
    #print(val)
    if (grepl("^Qualifier \\(", val)) {
      qualifier_transition_name <-
        val %>% stringr::str_replace("Qualifier", "QUAL") %>% stringr::str_replace("\\(", "") %>% stringr::str_replace("\\)", "\\]") %>% stringr::str_replace(" \\-\\> ", " \\> ")
      datWide[1, c] = paste(quantifier_name, " [", qualifier_transition_name, sep = "")
      col_name = datWide[1, c]
    } else if ((!is.na(val)) && nchar(val, keepNA = TRUE) > 0) {
      col_name = val
      quantifier_name = val
    } else if (stringr::str_detect(val, "Method")) {
      datWide[1, c] = datWide[1, c + 1]
      col_name = datWide[1, c]
    }
    else {
      datWide[1, c] = col_name
    }
  }
  datWide[1, ] <- as.list(stringr::str_squish(unlist(datWide[1, ])))

  # Concatenate rows containing parameters + transitions to the form parameter.sample and parameter.transition headers. This will later be converted with reshape()
  colnames(datWide) <- paste(datWide[2, ], datWide[1, ], sep = "\t")
  datWide <- datWide[-1:-2, ]

  # Rename some known column header names and remove columns that are not needed or not informative.
  old = c(
    "Data File\tSample",
    "Name\tSample",
    "AcqDate-Time\tSample",
    "Type\tSample",
    "Pos\tSample",
    "Vol\tSample"
  )
  new = c(
    "DataFileName",
    "SampleName",
    "AcqTimeStamp",
    "SampleType",
    "VialPosition",
    "InjVolume"
  )
  existing <- match(old, names(datWide))
  names(datWide)[na.omit(existing)] <- new[which(!is.na(existing))]

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
      AcqTimeStamp = lubridate::mdy_hm(.data$AcqTimeStamp),
      dplyr::across(.cols = dplyr::any_of(c("SampleName")),
                    .fns = stringr::str_squish)
    )

  # obtain list with column names of the columns that define transition values (e.g. "RT Cer d16:0/18:0"). Delimuter is currently tab (\t)
  param_transition_names <-
    colnames(datWide[, -1:-tail(grep("\\\t", colnames(datWide), invert =  TRUE), 1)])


  # Obtain long table of all param-transition combinations, split param and compund name and then spread values of different param as columns
  datLong <- datWide |>
    tidyr::pivot_longer(cols=param_transition_names, names_pattern = "(.*)\t(.*)$", names_to = c("Param", "Feature")) |>
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

  datLong <- datLong %>%
    dplyr::rename(Intensity = .data$Area,
                  DataName = .data$SampleName,
                  PrecursorMZ = .data$`Precursor Ion`,
                  ProductMZ = .data$`Product Ion`,
                  CollisionEnergy = .data$`Collision Energy`,
                  IonPolarity = .data$`Ion Polarity`)
  if(!silent) {
    cat("Imported ", length(unique(datLong$DataFileName)), "samples with ", fill = FALSE)
    cat(length(unique(datLong$Feature)), "transitions ", fill = FALSE)
  }

  datLong
}
