# RTL (Verilog) fully piplined implementation of integer FFT decimation-in-time algorithm
Fully piplined and axi-stream compatible implementation of integer FFT decimation-in-time algorithms in RTL (Verilog) with detailed explanation and math modelling (Python). The implementation is targeted at 7 Xilinx FPGA families (incl. Zynq 7000, Zynq MP) and tested on Zynq MPSoC. 
## Structure of the repository
Folders and files of the repository
```
+-- pptx\   
|   |
|   +-- fft_impl.pptx       -- presentation
|
+-- py_scripts\   
|   |
|   +-- data_generator.py   -- data generator for RTL simulation and fft_model.py
|   +-- fft_model.py        -- math model of FFT DIT algorithm for int16
|
+-- rtl\                    -- Verilog modules
|   |
|   +-- cmpx_adder.v
|   +-- cmpx_mult.v
|   +-- cmpx_reg.v
|   +-- FFT_DIT.v
|   +-- BRAM36E2.v
|   +-- f_stage.v
|   +-- l_stage.v
|   +-- n_stage.v
|   +-- radix2_bfly.v
|
+-- rtl_tb\              -- Testbench module
|   |
|   +-- FFT_DIT_tb.v
|
+-- tcl/                 -- TCL script for building Vivado project
|   |
|   +-- build.tcl
|
+-- Makefile             -- Makefile for building Vivado project
|
+-- README.md
|
+-- LICENSE

```

## Structure of the project
The project consists of set Verilog modules that describe Hardware Design of the FFT implementation based on DSP48E and BRAM36E2 blocks.

```
FFT_DIT.v - top module
|
+-- f_stage.v - the first stage of FFT algorithm
|   |
|   +--  radix2_bfly.v- radix 2 butterfly with decimation in time
|   |    |
|   |    +--  cmpx_mult.v  - complex multiplyer
|   |    +--  cmpx_adder.v - complex add / sub
|   |    +--  cmpx_reg.v   - complex register
|   |
|   +--  BRAM36E2.v
|
+-- n_stage.v - nth stage of FFT algorithm
|   |
|   +--  radix2_bfly.v- radix 2 butterfly with decimation in time
|   |    |
|   |    +--  cmpx_mult.v  - complex multiplyer
|   |    +--  cmpx_adder.v - complex add / sub
|   |    +--  cmpx_reg.v   - complex register
|   |
|   +--  BRAM36E2.v
|
+-- l_stage.v - the last stage of FFT algorithm for bit reversal or normal order
	|
	+--  BRAM36E2.v -- BRAM block instance
	|
	+--  BRAM36E2.v
```

## Top module parameters
| PARAMETER  | Meaning |
| ------------- |:-------------:|
| D_WIDTH      | data width for sample (16, 32, 64)     |
| N_RADIX      | amount of points (2**N_RADIX)     |
| D_ORDER      | output sorting ("NORMAL", "RESERVAL")     |