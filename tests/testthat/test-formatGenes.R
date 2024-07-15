# test script for formatGenes.R - testcases are NOT comprehensive!

test_that("formatGenes works", {
  data(apex_example.data)
  formatGenes_TenX <- formatGenes(immapex_example.data[["TenX"]],
                                    region = "v",
                                    technology = "TenX") 
  
  formatGenes_AIRR <- formatGenes(immapex_example.data[["AIRR"]],
                                    region = "v",
                                    technology = "AIRR", 
                                    simplify.format = FALSE) 
  
  formatGenes_OS <- formatGenes(immapex_example.data[["Omniscope"]],
                                  region = c("v", "j"),
                                  technology = "Omniscope") 

  formatGenes_Adaptive <- formatGenes(immapex_example.data[["Adaptive"]],
                                        region = "v",
                                        technology = "Adaptive") 
  
  expect_equal(
    formatGenes_Adaptive,
    getdata("formatGenes", "formatGenes_Adaptive")
  )
  
  
  expect_equal(
    formatGenes_TenX,
    getdata("formatGenes", "formatGenes_TenX")
  )
  
  
  expect_equal(
    formatGenes_OS,
    getdata("formatGenes", "formatGenes_OS")
  )
  
  
  expect_equal(
    formatGenes_AIRR,
    getdata("formatGenes", "formatGenes_AIRR")
  )
})

