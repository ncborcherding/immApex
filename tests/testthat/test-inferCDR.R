# test script for inferCDR.R - testcases are NOT comprehensive!

test_that("inferCDR works", {
  
  data(immapex_example.data)
  reference <- getdata("getIMGT", "getIMGT_TRBV_human_inframe_aa")
  
  TenX_formatted <- formatGenes(immapex_example.data[["TenX"]],
                                region = "v",
                                technology = "TenX")
  
  TenX_formatted <- inferCDR(TenX_formatted,
                             chain = "TRB", 
                             reference = reference,
                             technology = "TenX", 
                             sequence.type = "aa",
                             sequences = c("CDR1", "CDR2"))[,c("CDR1_IMGT","CDR2_IMGT")]
  
  expect_equal(
    TenX_formatted,
    getdata("inferCDR", "inferCDR_Tenx")
  )
  
  AIRR_formatted <- formatGenes(immapex_example.data[["AIRR"]],
                                region = "v",
                                technology = "AIRR")
  
  AIRR_formatted <- inferCDR(AIRR_formatted,
                              chain = "TRB", 
                              reference = reference,
                              technology = "TenX", 
                              sequence.type = "aa",
                              sequences = c("CDR1"))[,"CDR1_IMGT"]
  
  expect_equal(
    AIRR_formatted,
    getdata("inferCDR", "inferCDR_AIRR")
  )
  
})

