# test script for formatGenes.R - testcases are NOT comprehensive!

test_that("formatGenes works", {
  data(immapex_example.data)
  formatGenes_TenX <- formatGenes(immapex_example.data[["TenX"]],
                                    region = "v",
                                    technology = "TenX")[,c("v_IMGT", "v_IMGT.check")]
  
  formatGenes_AIRR <- formatGenes(immapex_example.data[["AIRR"]],
                                    region = "v",
                                    technology = "AIRR", 
                                    simplify.format = FALSE)[,c("v_IMGT", "v_IMGT.check")]
  

  formatGenes_Adaptive <- formatGenes(immapex_example.data[["Adaptive"]],
                                        region = "v",
                                        technology = "Adaptive")[,c("v_IMGT", "v_IMGT.check")]
  
  expect_equal(
    formatGenes_Adaptive,
    getdata("formatGenes", "formatGenes_Adaptive")
  )
  
  
  expect_equal(
    formatGenes_TenX,
    getdata("formatGenes", "formatGenes_TenX")
  )
  
  expect_equal(
    formatGenes_AIRR,
    getdata("formatGenes", "formatGenes_AIRR")
  )
})

