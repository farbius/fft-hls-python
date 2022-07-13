############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
############################################################
open_project -reset fft_accel
set_top FFT_TOP
add_files coef_init.h
add_files fft_accel.cpp
add_files fft_accel.h
add_files -tb fft_accel_tb.cpp
add_files -tb testbenches.h


open_solution -reset "solution1"
set_part {xczu2cg-sfvc784-1-e}
create_clock -period 10 -name default
#source "./fft_accel/solution1/directives.tcl"
csim_design

exit

csynth_design
cosim_design
export_design -format ip_catalog
