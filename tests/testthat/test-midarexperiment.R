testthat::test_that("import_MSOrganizerXLM: Template with all information", {
  mexp <- SLINGtools::MidarExperiment()
  mexp <- SLINGtools::loadMasshunterCSV(mexp, "21_Test_MH.csv")
  mexp <- SLINGtools::loadMSOrganizerXLM(mexp, "20_MSTemplate_Creator_forTest.xlsm")

  mexp <- SLINGtools::normalizeByISTD(mexp)
  mexp <- SLINGtools::quantitateByISTD(mexp)
  mexp <- SLINGtools::calcQC(mexp)

  dd <- readRDS("21_MidarExperiment_1.rds")
  expect_equal(dplyr::all_equal(mexp@d_QC, dd@d_QC), TRUE) &
  expect_equal(dplyr::all_equal(mexp@dataset, dd@dataset), TRUE)
})

