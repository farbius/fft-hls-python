# High-Level Synthesis fully piplined implementation of integer FFT decimation-in-time algorithm
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
|   +-- fft-hls-python.md   -- main documentation 
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
|   +-- coef_init.h
|
+-- tcl\              		-- tcl script for compiling Vivado HLS
|   |
|   +-- run_hls.tcl
|
+-- Makefile                -- Makefile for building Vivado HLS project
|
+-- README.md
|
+-- LICENSE
```

How to work with the repository

1. Create synthetic signal with the <i>signal_generator.py</i> script (explained [here](./doc/fft-hls-python.md#mathematical-modelling))

```sh
python3 py_scripts/signal_generator.py 1024 3 40
```

2. In the repository root directory run <i>make</i> command for building and launching Vivado HLS project

```sh
make
```

3. Run <i>fft_model.py</i> script for cheking Vivado HLS simulation result

```sh
python3 py_scripts/fft_model.py
```
