# test script for infer.cdr.R - testcases are NOT comprehensive!

test_that("infer.cdr works", {
  
  data(apex_example.data)
  reference <- getdata("get.IMGT", "get.IMGT_TRBV_human_inframe_aa")
  
  TenX_formatted <- format.genes(apex_example.data[["TenX"]],
                                 region = "v",
                                 technology = "TenX")
  
  TenX_formatted <- infer.cdr(TenX_formatted,
                              chain = "TRB", 
                              reference = reference,
                              technology = "TenX", 
                              sequence.type = "aa",
                              sequences = c("CDR1", "CDR2"))
  
  expect_equal(
    TenX_formatted,
    getdata("infer.cdr", "infer.cdr_Tenx")
  )
  
  AIRR_formatted <- format.genes(apex_example.data[["AIRR"]],
                                 region = "v",
                                 technology = "AIRR")
  
  AIRR_formatted <- infer.cdr(AIRR_formatted,
                              chain = "TRB", 
                              reference = reference,
                              technology = "TenX", 
                              sequence.type = "aa",
                              sequences = c("CDR1"))
  
  expect_equal(
    AIRR_formatted,
    getdata("infer.cdr", "infer.cdr_AIRR")
  )
  
})

