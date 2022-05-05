



#' Read and convert an Agilent MassHunter Quant CSV result file
#'
#' @param rawFileName
#'
#' @return A tibble in the long format
#' @export
#'
#' @examples
#' library(SLINGtools)
#'
#' data_file_path <- system.file("extdata", "Testdata_Lipidomics_MHQuant_Detailed.csv", package = "SLINGtools")
#' d <- read_MassHunterCSV(data_file_path)
#' d
#'
read_MassHunterCSV <- function(rawFileName) {

    cat(paste0("Importing ", rawFileName, ": "), fill = FALSE)

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


  # Fill in compound name in empty columns (different parameters of the same compound)
  # Dealing with exported "Qualifier" results: Adds the corresponding quantifier name (the first column of the group) with the Tag "QUAL" and the transition in front e.g. "Sph d18:0 [QUAL 302.3>266.2]"
  col_name = "Sample"
  isFirstQualifier = TRUE
  for (c in 1:ncol(datWide)) {
    val = datWide[1, c]
    #print(val)
    if (grepl("^Qualifier \\(", val)) {
      qualifier_transition_name <-
        val %>% stringr::str_replace(., "Qualifier", "QUAL") %>% stringr::str_replace(., "\\(", "") %>% stringr::str_replace(., "\\)", "\\]") %>% stringr::str_replace(., " \\-\\> ", " \\> ")
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


  # Concatenate rows containing parameters + compounds to the form parameter.sample and parameter.compound headers. This will later be converted with reshape()
  datWide <- datWide %>% {
    setNames(., paste(.[2, ], .[1, ], sep = "\t"))
  }
  datWide <-
    datWide[-1:-2, ]  # remove lines with header and first two empty columns (Sample and empty with only !)



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


  # obtain list with column names of the columns that define compound values (e.g. "RT Cer d16:0/18:0"). Delimuter is currently tab (\t)
  param_compound_names <-
    colnames(datWide[, -1:-tail(grep("\\\t", colnames(datWide), invert =  TRUE), 1)])

  # Obtain long table of all param-compound combinations, split param and compund name and then spread values of different param as columns
  datLong <- datWide %>%
    tidyr::gather(key = "ParamCompound",
           value = "Value",
           all_of(param_compound_names)) %>%
    tidyr::separate(col = ParamCompound, into = c("Param", "Compound"), "\t") %>%
    tidyr::spread(Param, Value)

  # Convert types of knows parameters and fields in the data set
  # ------------------------------------------------------------
  datLong <- datLong %>%
    dplyr::mutate(Compound = trimws(Compound)) %>%
    dplyr::mutate_if(is.character, stringr::str_squish) %>%
    dplyr::mutate_at(.vars = dplyr::vars(
      dplyr::matches(
        "RT|Area|Height|FWHM|Width|SN|IntStart|IntEnd|Symmetry|InjVolume"
      )
    ),
    .funs = list( ~ as.numeric(.))) %>%
    dplyr::mutate_at(.vars = dplyr::vars(
      dplyr::matches(
        "Compound|QuantWarning|DataFileName|SampleName|SampleType|VialPosition|Method"
      )
    ),
    .funs = list( ~ stringr::str_squish(.))) %>%
    dplyr::mutate_at(.vars = dplyr::vars(matches("AcqTimeStamp")), .funs = list( ~ lubridate::dmy_hm(.)))
  datLong <- datLong %>% dplyr::ungroup() %>% dplyr::rename(Intensity = Area,
                                              DataName = SampleName)
  cat(length(unique(datLong$Compound)), "transitions; ", fill = FALSE)
  cat(length(unique(datLong$DataFileName)), "samples", fill = TRUE)
  return(datLong)
}
