set ProjectName  "fft_hw"
set PartDev      "xczu2cg-sfvc784-1-i" 

set TclPath 	[file dirname [file normalize [info script]]]
set TclPath 	[file dirname $TclPath]
set ProjectPath  $TclPath/build


open_project fft_accel
set_top FFT_TOP
add_files ../hls_src/fft_accel.h
add_files ../hls_src/fft_accel.cpp
add_files -tb ../hls_src/testbenches.h -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
add_files -tb ../hls_src/fft_accel_tb.cpp -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
open_solution "solution1" -flow_target vivado
set_part {xczu4ev-sfvc784-1-e}
create_clock -period 10 -name default
