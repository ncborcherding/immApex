# test script for getIMGT.R - testcases are NOT comprehensive!

test_that("getIMGT works", {

  #Default Test
  TRBV_human_inframe_aa <- getIMGT(species = "human",
                                    chain = "TRB",
                                    frame = "inframe",
                                    region = "v",
                                    sequence.type = "aa") 
  
  expect_equal(
    TRBV_human_inframe_aa,
    getdata("getIMGT", "getIMGT_TRBV_human_inframe_aa")
  )
  
  #Test Different Region and Species
  TRAJ_mouse_inframe_aa <- getIMGT(species = "mouse",
                                    chain = "TRB",
                                    frame = "inframe",
                                    region = "j",
                                    sequence.type = "aa") 
  
  expect_equal(
    TRAJ_mouse_inframe_aa,
    getdata("getIMGT", "getIMGT_TRAJ_mouse_inframe_aa")
  )
  
  #Test All Sequence Pull
  IGHV_rat_all_nt <- getIMGT(species = "rat",
                              chain = "IGH",
                              frame = "all",
                              region = "v",
                              sequence.type = "nt") 
  
  expect_equal(
    IGHV_rat_all_nt,
    getdata("getIMGT", "getIMGT_IGHV_rat_all_nt")
  )
  
  #Test IMGT Gap Sequence Pull
  TRBV_rabbit_gap_aa <- getIMGT(species = "rabbit",
                              chain = "TRB",
                              frame = "inframe+gap",
                              region = "v",
                              sequence.type = "aa") 
  
  expect_equal(
    TRBV_rabbit_gap_aa,
    getdata("getIMGT", "getIMGT_TRBV_rabbit_gap_aa")
  )
  
  
  TRAJ_pig_inframe_aa <- getIMGT(species = "pig",
                                 chain = "TRA",
                                 frame = "inframe",
                                 region = "v",
                                 sequence.type = "aa") 
  
  expect_equal(
    TRAJ_pig_inframe_aa,
    getdata("getIMGT", "getIMGT_TRAJ_pig_inframe_aa")
  )
  
  
  expect_equal(
    .parseSpecies("ferret"),
    "Mustela putorius furo"
  )
  
  expect_equal(
    .parseSpecies("Ferret"),
    "Mustela putorius furo"
  )
  
  expect_equal(
    .parseSpecies("rhesus monkey"),
    "Macaca mulatta"
  )
  
  expect_equal(
    .parseSpecies("Rhesus monkey"),
    "Macaca mulatta"
  )
  
  expect_equal(
    .parseSpecies("rhesus Monkey"),
    "Macaca mulatta"
  )
  
})

