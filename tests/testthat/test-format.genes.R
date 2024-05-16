# test script for format.genes.R - testcases are NOT comprehensive!

test_that("format.genes works", {

  format.genes_TenX <- format.genes(Apex_example.data[["TenX"]],
                                    region = "v",
                                    technology = "TenX") 
  
  format.genes_AIRR <- format.genes(Apex_example.data[["AIRR"]],
                                    region = "v",
                                    technology = "AIRR", 
                                    simplify.format = FALSE) 
  
  format.genes_OS <- format.genes(Apex_example.data[["Omniscope"]],
                                  region = c("v", "j"),
                                  technology = "Omniscope") 

  format.genes_Adaptive <- format.genes(Apex_example.data[["Adaptive"]],
                                        region = "v",
                                        technology = "Adaptive") 
  
  expect_equal(
    format.genes_Adaptive,
    getdata("format.genes", "format.genes_Adaptive")
  )
  
  
  expect_equal(
    format.genes_TenX,
    getdata("format.genes", "format.genes_TenX")
  )
  
  
  expect_equal(
    format.genes_OS,
    getdata("format.genes", "format.genes_OS")
  )
  
  
  expect_equal(
    format.genes_AIRR,
    getdata("format.genes", "format.genes_AIRR")
  )
})

