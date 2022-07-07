# RTL (Verilog) fully piplined implementation of integer FFT decimation-in-time algorithm
Fully piplined and axi-stream compatible implementation of integer FFT decimation-in-time algorithms in Vivado HLS with detailed explanation and math modelling (Python).


## Structure of the repository
Folders and files of the repository
```
+-- doc\   
|   |
|   +-- fft-hls-python.md       -- main doc 
|   +-- 1.md       				-- main doc 
|   +-- 2.md       				-- main doc 
|   +-- 3.md       				-- main doc
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
