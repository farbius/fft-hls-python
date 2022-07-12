# RTL (Verilog) fully piplined implementation of integer FFT decimation-in-time algorithm
Fully piplined and axi-stream compatible implementation of integer FFT decimation-in-time algorithms in Vivado HLS with detailed explanation and math modelling (Python).

Explanation

1. [Theory](./doc/fft-hls-python.md#theory)
2. [Mathematical Modelling](./doc/fft-hls-python.md#mathematical-modelling)
3. [High Level Synthesis implementation with C/C++](./doc/fft-hls-python.md#high-level-synthesis-implementation)
4. [References](./doc/fft-hls-python.md#references)

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
|   +-- signal_generator.py -- complex data generator for RTL simulation and fft_model.py
|   +-- fft_model.py        -- math model of FFT DIT algorithm for int16
|
+-- hls_src\                -- Vivado HLS sources
|   |
|   +-- fft_accel.cpp
|   +-- fft_accel.h
|   +-- fft_accel_tb.cpp
|   +-- testbenches.h
|
+-- sim_files\              -- input / output files
|   |
|   +-- cmpx_hls.txt
|   +-- nonscaled_re.txt
|   +-- nonscaled_im.txt
|   +-- scaled_re.txt
|   +-- scaled_im.txt
|
+-- Makefile                -- Makefile for building Vivado HLS project
|
+-- README.md
|
+-- LICENSE
```
