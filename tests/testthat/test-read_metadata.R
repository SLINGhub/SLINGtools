testthat::test_that("import_MSOrganizerXLM: Template with all information", {
    d <- SLINGtools::import_MSOrganizerXLM("20_MSTemplate_Creator_forTest.xlsm")
    dd <- readRDS("20_MSTemplate_Creator_forTest.rds")
    expect_identical(d, dd)
})

