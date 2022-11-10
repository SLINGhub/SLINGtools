
# update all test reference data... I/you know what you are doing...


d <- SLINGtools::read_MassHunterCSV("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv")
saveRDS(object = d,  version = 2, file = "1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.rds")


d <- SLINGtools::read_MassHunterCSV("2_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM.csv")
saveRDS(object = d,  version = 2, file = "2_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM.rds")


d <- SLINGtools::read_MassHunterCSV("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv")
saveRDS(object = d,  version = 2, file = "3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.rds")

d <- SLINGtools::read_MassHunterCSV("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv")
saveRDS(object = d,  version = 2, file = "4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.rds")

d <- SLINGtools::read_MassHunterCSV("5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.csv")
saveRDS(object = d,  version = 2, file = "5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.rds")

d <- SLINGtools::read_MassHunterCSV("6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.csv")
saveRDS(object = d,  version = 2, file = "6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.rds")

d <- SLINGtools::read_MassHunterCSV("7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.csv")
saveRDS(object = d,  version = 2, file = "7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.rds")

d <- SLINGtools::read_MassHunterCSV("8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.csv")
saveRDS(object = d,  version = 2, file = "8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.rds")

d <- SLINGtools::read_MassHunterCSV("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.csv")
saveRDS(object = d,  version = 2, file = "9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.rds")

d <- SLINGtools::read_MassHunterCSV("10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.csv")
saveRDS(object = d,  version = 2, file = "10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.rds")

d <- SLINGtools::read_MassHunterCSV("11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.csv")
saveRDS(object = d,  version = 2, file = "11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.rds")

d <- SLINGtools::read_MassHunterCSV("15_Testdata_MHQuant_Corrupt_ExtraTopLine.csv")
saveRDS(object = d,  version = 2, file = "15_Testdata_MHQuant_Corrupt_ExtraTopLine.rds")

d <- SLINGtools::read_MassHunterCSV_wide("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv", field = "Area")
saveRDS(object = d,  version = 2, file = "3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults-Area-flat.rds")

d <- SLINGtools::read_MassHunterCSV("17_Testdata_Lipidomics_GermanSystem.csv")
saveRDS(object = d,  version = 2, file = "17_Testdata_Lipidomics_GermanSystem.rds")


d <- SLINGtools::read_MassHunterCSV("18_Testdata_MHQuant_DefaultSampleInfo_AreaOnly_AMPM.csv")
saveRDS(object = d,  version = 2, file = "18_Testdata_MHQuant_DefaultSampleInfo_AreaOnly_AMPM.rds")

