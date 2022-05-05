testthat::test_that("Import Detailed MH CSV", {

    d <- SLINGtools::read_MassHunterCSV("Testdata_Lipidomics_MHQuant_Detailed.csv")
    dd <- readRDS("Testdata_Lipidomics_MHQuant_Detailed.rds")
    expect_identical(d, dd)
})
