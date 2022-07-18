
# <div align="center">Fast Fourier Transform mathematical modelling and HLS implementation</div>


<div align="right"> <i>PhD, Senior R&D Engineer </i></div>
<div align="right"> <i>Aleksei Rostov </i></div>
<div align="right"> <i>Munich, 2022</i> </div>

## Agenda
1. [Theory](#theory)
2. [Mathematical Modelling](#mathematical-modelling)
3. [High Level Synthesis implementation with C/C++](#high-level-synthesis-implementation)
4. [References](#references)




## Theory
Discrete Fourier Transform (DFT) of a finite-length sequence of length *N* is

$$X[k]  =\sum_{n=0}^{N-1}x[n] W_{N}^{kn}, \qquad k = 0,1,...,N-1.\qquad (1.1)$$

where $W_{N}^{kn} = e^{-j(2\pi kn/N)}$.

<p align="justify">
Direct computation of $X[k]$ requires a total of $N^2$ complex multiplications and $N(N-1)$ complex additions.
Fast Fourier Transform (FFT)  is <i>exactly</i> the same DFT with optimization by reducing number of computations.
All optimizations for improving the efficiency of the computation are based on the symmetry and periodicity properties of   $W_{N}^{kn}$ [1], specifically,

$$W_{N}^{k[N - n]} =W_{N}^{-kn} = (W_{N}^{kn})^* \qquad (symmetry ) \qquad (1.2)$$


$$W_{N}^{kn} =W_{N}^{k(n+N)} = W_{N}^{(k+N)n} \qquad (periodicity) \qquad (1.3)$$


For explanation let's consider direct calculation of two samples of  $X[k]$ from Eq. (1.1) for $N=8$

$$X[2]  =x[0] W_{8}^{0} + x[1] W_{8}^{2} + x[2] W_{8}^{4}+x[3] W_{8}^{6} + x[4] W_{8}^{8} + x[5] W_{8}^{10}+ x[6] W_{8}^{12} + x[7] W_{8}^{14} \qquad$$

$$X[3]  =x[0] W_{8}^{0} + x[1] W_{8}^{3} + x[2] W_{8}^{6}+x[3] W_{8}^{9} + x[4] W_{8}^{12} + x[5] W_{8}^{15}+ x[6] W_{8}^{18} + x[7] W_{8}^{21} \qquad$$

By using the periodicity property of   $W_{N}^{kn}$ and the fact that  $W_{N}^{N/2} =  e^{-j(2\pi/N)N/2}=-1$, we obtain

$$X[2]  =x[0] W_{8}^{0} + x[1] W_{8}^{2} - x[2] W_{8}^{0}-x[3] W_{8}^{2} + x[4] W_{8}^{0} + x[5] W_{8}^{2}- x[6] W_{8}^{0} - x[7] W_{8}^{2} \qquad$$

$$X[3]  =x[0] W_{8}^{0} + x[1] W_{8}^{3} - x[2] W_{8}^{2}+x[3] W_{8}^{1} - x[4] W_{8}^{0} - x[5] W_{8}^{3}+ x[6] W_{8}^{2} - x[7] W_{8}^{1} \qquad$$

<p align="justify">
By reducing number of $W_{N}^{kn}$ we may group corresponding pairs of $x[n]$  samples with equal coefficients, that allow to decrease number of multiplications. Also coefficients transformation can be explained by Fig. (1.1)

<p align="center">
  <img src="https://github.com/farbius/fft-hls-python/blob/main/doc/images/w_phasor.png" alt="w_phasor"/>
</p>

<div align="center">
<b>Figure 1.1 </b> Transformation of coefficients
</div>
<br/>

As seen from Fig. (1.1) , it is enough to compute and store only $N/2$ of coefficients. This is the main advantage of periodocity property of $W_{N}^{kn}$. By applying both the symmetry and periodicity property  for coefficients and grouping input data the following dataflow for 8-point DFT computation can be implemented

<p align="center">
  <img src="https://github.com/farbius/fft-hls-python/blob/main/doc/images/data_flow.png" alt="data_flow"/>
</p>

<div align="center">
<b>Figure 1.2 </b> Flow graph of decomposition and computation of 8-point DTF
</div>
<br/> 

<p align="justify">
This implementation requires number of multiplications and additions is equal  to $N\log_2N$ , what is significantly less then for direct implementation of DFT $N^2$. The reduction extremelly grows up with the increasing number of DFT points. (Here I'm focused only on decimation-in-time algorithm of FFT, but of course there is amount of other algorithms for implementation FFT).
FFT implementation consist of $\log_2N$ stages and $N/2$ butterflies for every stage. The butterfly is a unit of FFT that implemets 2-point DFT with  one  $W_{N}^{kn}$ coefficient  Fig. (1.3)

<p align="center">
  <img src="https://github.com/farbius/fft-hls-python/blob/main/doc/images/butterfly.png" alt="butterfly"/>
</p>

<div align="center">
<b>Figure 1.3 </b> Flow graph of basic butterfly computation
</div>
<br/> 

Associated equations for a butterfly from Fig. (1.3)

$$X_m[p] =X_{m-1}[p] + W_N^rX_{m-1}[q]  \qquad (1.4)$$
$$X_m[q] =X_{m-1}[p] - W_N^rX_{m-1}[q]  \qquad (1.5)$$

<p align="justify">
The batterfly requires only one complex multiplication $W_N^rX_{m-1}[q]$ and $N\log_2N$ multiplication for computing FFT overall. For example, 8-point FFT requires 24  complex multipliers as can be seen from Fig. (1.2).
Apart of butterfly calculation it is necessary to reoder input data for the first stage. 
Reordering involves <i>bit-reversal</i> algorithm, when every bit of input data index binary form is reversed. For example, indexes for  $N = 8$ will be reorder like

<p align="justify">

| Normal Index| Binary form of normal index |Binary form of reserved index|Reserved index|
| :-: | :-:| :-:| :-:|
| 0 | 000 | 000 | 0 |
| 1 | 001 | 100 | 4 |
| 2 | 010 | 010 | 2 |
| 3 | 011 | 110 | 6 |
| 4 | 100 | 001 | 1 |
| 5 | 101 | 101 | 5 |
| 6 | 110 | 011 | 3 |
| 7 | 111 | 111 | 7 |



Summarize the information above, FFT implementation requires (decimation-in-time Cooley–Tukey FFT algorithm):
1. butterfly for complex multiplication on every stage
1. set of complex coefficients
1. input signal reodering



<br/>

## Mathematical Modelling


<p align="justify">
Mathematical modelling of FFT allows to figure out <i>hardware design</i> of the algorithm and proof outputs during HLS Co-Simulation. 
<br>
For mathematical modelling two Python scripts were created: <i>signal_generator.py</i> and <i>fft_model.py</i>.  First script generates a synthetic complex signal for HLS testbench and header <i>coef_init.h</i> file for initialization of HLS FFT implementation. The script may be launching with, for example,the following arguments

```sh
python3 signal_generator.py 1024 3 40
```
what means - the synthetic signal consist of 1024 samples and is composed from 3 signals with Signal-To-Noise ratios 40 dB, 38 dB, 36 dB. 
Code for parsing the arguments is 

```sh
parser.add_argument('Npoints' ,default=None, type=int)
parser.add_argument('Nsignals',default=None, type=int)
parser.add_argument('SNR_dB'  ,default=None, type=int)
```
<p align="justify">
where <i>Npoints</i>  is a number of FFT points, 
<i>Nsignals</i> - is an amount of harmonics with random frequencies in the synthetic signal and <i>SNR_dB</i> - is a max value of Signal-To-Noise ratio of the first harmonic in the synthetic signal, SNR of other harmonics will be decreased on 2 dB for every one.
Mathematical model of the complex synthetic signal is

$$x_k[n] = 10^{(SNR_{dB}-2k)/20}e^{-j2\pi f_k n/Npoints} + n_k[n] \qquad (2.1)$$
<p align="justify">
where $n = 0 ... Npoints-1$ - signal's sample, $k = 0 ... Nsignals-1$ - harmonic's number, $n_k$ - random noise for every harmonic.
Frequency meaning $f_k$ of harmonic  is choosen randomly and a power is decreased for every harmonic on 2 dB.
<br>
After calculation in Eq. (2.1) the amplitude of synthetic sum signal is scaled to <i>-1 ... +1</i> and <i>-32768 ... +32768</i> and samples are saved to files.
<br>
The script will generate several files:

```sh
nonscaled_re.txt  
nonscaled_im.txt
scaled_re.txt
scaled_im.txt
coef_init.h
```
<p align="justify">
Header <i>coef_init.h</i> consist of FFT parameters and scaled to  <i>-32768 ... +32768</i> coefficients. The header file is used by HLS FFT implementation. 
Files <i>nonscaled_.txt</i> are complex float point input for Numpy FFT implementation, that will be used for comparison with HLS FFT implementation. 
Files <i>scaled_.txt</i> are scaled to 16-bit signed register complex input for HLS FFT.
<br>
The <i>fft_model.py</i> python script consists of scaled to 16-bits signed register implementation of FFT and FFT from Numpy package. The script reads out <i>scaled_.txt</i> and <i>nonscaled_.txt</i> data, output <i>cmpx_hls.txt</i> from HLS Co-simulation and plots results of three FFT. 
<br>
Result for <i>python3 signal_generator.py 1024 3 40</i> is depicted on Fig. (2.1)


<p align="center">
  <img src="https://github.com/farbius/fft-hls-python/blob/main/doc/images/fft_snr.png" alt="butterfly"/>
</p>

<div align="center">
<b>Figure 2.1 </b> Results of Numpy, Python and HLS FFT implementation
</div>
<br/> 
<p align="justify">
First plot is a power of FFT output implemented in Numpy (green), Python (red) and HLS Co-simulation (blue). Since Python and HLS implementation is scaled to 16-bit signed register, thir nois floor is restricted by Root Means Square (RMS) value, that depends on length of a register. In case of 16-bit signed register, minimum level of power spectrum can be 

$$SNR_{RMS} = 6.02 * 15 + 1.76 = 92.06 \  dB$$

<p align="justify">
FFT itself without scaling, expands power spectrum range on value 

$$SNR_{FFT} = 10*\log_{10}(Npoints)  \  dB$$
<p align="justify">
Noise floor reduction in FFT is caused by narrow-bandness  of FFT itself [2]. For demonstration the noise floor reduction, let's consider two noise 1024-points  signals with amplitude 0 dB and -20 dB in the following Python script

```sh
import numpy as np
import matplotlib.pyplot as plt

Npoints = 1024
SNR_1 = 0
SNR_2 = -20

n_1 = 10**(SNR_1/20)*np.random.randn(Npoints)
n_2 = 10**(SNR_2/20)*np.random.randn(Npoints)

# normalized Numpy FFT 
nf_1 = np.abs(np.fft.fft(n_1)/Npoints)
nf_2 = np.abs(np.fft.fft(n_2)/Npoints)

plt.figure()
plt.plot(20*np.log10(nf_1), '.-r', label='nf_1 SNR = 0 dB')
plt.plot(20*np.log10(nf_2), '.-b', label='nf_2 SNR = -20 dB')
plt.title("FFT {} points, FFT noise floor is {:3.2f} dB".format(Npoints, 10*np.log10(Npoints)), fontweight="bold", fontsize=14)
plt.legend(loc='upper right')
plt.grid()
plt.xlabel('bin')
plt.ylabel('power, dB')
plt.show()
```

<p align="center">
  <img src="https://github.com/farbius/fft-hls-python/blob/main/doc/images/fft_floor_snr.png" alt="butterfly"/>
</p>

<div align="center">
<b>Figure 2.2 </b> FFT noise floor reduction
</div>
<br/> 
<p align="justify">
As can be seen from Fig. (2.2)  noise floor was moved on 30 dB since 1024-point FFT was applied (30 dB).
<br>
<p align="justify">
Second plot in Fig. (2.1) is a result of comparison Python FFT implementation and HLS Co-simulation output. The small error is caused by rounding operation during FFT computation. 




## High Level Synthesis Implementation





## References

1. Alan V Oppenheim, Ronald W. Schafer, Discrete-Time Signal Processing, 3rd Edition, 2010
2. [W. Kester, Understand SINAD, ENOB, SNR, THD, THD + N, and SFDR so You Don't Get Lost in the Noise Floor, Analog Devices, 2008](https://www.analog.com/media/en/training-seminars/tutorials/MT-003.pdf)
3. source