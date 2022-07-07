# RTL (Verilog) fully piplined implementation of integer FFT decimation-in-time algorithm
Fully piplined and axi-stream compatible implementation of integer FFT decimation-in-time algorithms in Vivado HLS with detailed explanation and math modelling (Python).

Explanation can be found [here](https://github.com/farbius/fft-hls-python/tree/main/doc)
## Structure of the repository
Folders and files of the repository
```
+-- doc\   
|   |
|   +-- fft-hls-python.md       -- main documentation 
|	+-- images\
|
+-- py_scripts\   
|   |
|   +-- data_generator.py   -- data generator for RTL simulation and fft_model.py
|   +-- fft_model.py        -- math model of FFT DIT algorithm for int16
|
+-- hls_src\                -- Vivado HLS sources
|   |
|   +-- fft_accel.cpp
|   +-- fft_accel.h
|   +-- fft_accel_tb.cpp
|   +-- testbenches.cpp
|
+-- sim_files\              -- Files for simulation
|   |
|   +-- FFT_DIT_tb.v
|
+-- tcl/                    -- TCL script for building Vivado HLS project
|   |
|   +-- build.tcl
|
+-- Makefile                -- Makefile for building Vivado HLS project
|
+-- README.md
|
+-- LICENSE
```
