# RTL (Verilog) fully piplined implementation of integer FFT decimation-in-time algorithm
Fully piplined and axi-stream compatible implementation of integer FFT decimation-in-time algorithms in Vivado HLS with detailed explanation and math modelling (Python).

Explanation can be found [here](https://github.com/farbius/fft-hls-python/blob/main/doc/fft-hls-python.md)
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
