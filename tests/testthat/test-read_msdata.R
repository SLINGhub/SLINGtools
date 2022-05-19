testthat::test_that("read_MassHunterCSV: Detailed MH CSV", {

    d <- SLINGtools::read_MassHunterCSV("Testdata_Lipidomics_MHQuant_Detailed.csv")
    dd <- readRDS("Testdata_Lipidomics_MHQuant_Detailed.rds")
    expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: Detailed MH CSV but AcqDateTime missing", {

  d <- SLINGtools::read_MassHunterCSV("Testdata_Lipidomics_MHQuant_Short_Detailed_NoAcqDateTime.csv")
  dd <- readRDS("Testdata_Lipidomics_MHQuant_Short_Detailed_NoAcqDateTime.rds")
  expect_identical(d, dd)
})



testthat::test_that("read_peakareas_MassHunterCSV: Detailed MH CSV", {
  d <- SLINGtools::read_peakareas_MassHunterCSV("Testdata_Lipidomics_MHQuant_Detailed.csv")
  dd <- readRDS("Areas_Testdata_Lipidomics_MHQuant_Detailed.rds")
  expect_identical(d, dd)
})
