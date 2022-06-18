set ProjectName  "fft_hw"
set PartDev      "xczu2cg-sfvc784-1-i" 

set TclPath 	[file dirname [file normalize [info script]]]
set TclPath 	[file dirname $TclPath]
set ProjectPath  $TclPath/build

create_project $ProjectName $ProjectPath -part $PartDev

add_files -norecurse $TclPath/rtl/BRAM36E2.v
add_files -norecurse $TclPath/rtl/cmpx_reg.v
add_files -norecurse $TclPath/rtl/FFT_DIT.v
add_files -norecurse $TclPath/rtl/l_stage.v
add_files -norecurse $TclPath/rtl/cmpx_adder.v
add_files -norecurse $TclPath/rtl/n_stage.v
add_files -norecurse $TclPath/rtl/radix2_bfly.v
add_files -norecurse $TclPath/rtl/f_stage.v
add_files -norecurse $TclPath/rtl/cmpx_mult.v

set_property SOURCE_SET sources_1 [get_filesets sim_1]
add_files -fileset sim_1 -norecurse $TclPath/rtl_tb/FFT_DIT_tb.v

set_property top FFT_DIT_tb [get_filesets sim_1]
set_property top_lib xil_defaultlib [get_filesets sim_1]
update_compile_order -fileset sim_1
