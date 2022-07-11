
# <div align="center">Fast Fourier Transform mathematical modelling and HLS implementation</div>


<div align="right"> <i>PhD, Senior R&D Engineer </i></div>
<div align="right"> <i>Aleksei Rostov </i></div>
<div align="right"> <i>Munich, 2022</i> </div>

## Agenda
1. [Theory](#theory)
2. [Mathematical modelling and SNR explanation](#mathematical)
3. [High Level Synthesis implementation with C/C++](#high)
4. [References](#references)




## Theory
Discrete Fourier Transform (DFT) of a finite-length sequence of length *N* is

$$X[k]  =\sum_{n=0}^{N-1}x[n] W_{N}^{kn}, \qquad k = 0,1,...,N-1.\qquad (1.1)$$

where $W_{N}^{kn} = e^{-j(2\pi kn/N)}$.
Direct computation of $X[k]$ requires a total of $N^2$ complex multiplications and $N(N-1)$ complex additions.
Fast Fourier Transform (FFT)  is *exactly* the same DFT with optimization by reducing number of computations.
All optimizations for improving the efficiency of the computation are based on the symmetry and periodicity ptoperties of   $W_{N}^{kn}$, specifically,

$$W_{N}^{k[N - n]} =W_{N}^{-kn} = (W_{N}^{kn})^* \qquad (symmetry ) \qquad (1.2)$$
$$W_{N}^{kn} =W_{N}^{k(n+N)} = W_{N}^{(k+N)n} \qquad (periodicity) \qquad (1.3)$$

For explanation let's consider direct calculation of two samples of  $X[k]$ from Eq. (1.1) for $N=8$

$$X[2]  =x[0] W_{8}^{0} + x[1] W_{8}^{2} + x[2] W_{8}^{4}+x[3] W_{8}^{6} + x[4] W_{8}^{8} + x[5] W_{8}^{10}+ x[6] W_{8}^{12} + x[7] W_{8}^{14} \qquad$$
$$X[3]  =x[0] W_{8}^{0} + x[1] W_{8}^{3} + x[2] W_{8}^{6}+x[3] W_{8}^{9} + x[4] W_{8}^{12} + x[5] W_{8}^{15}+ x[6] W_{8}^{18} + x[7] W_{8}^{21} \qquad$$

By using the periodicity property of   $W_{N}^{kn}$ and the fact that  $W_{N}^{N/2} =  e^{-j(2\pi/N)N/2}=-1$, we obtain

$$X[2]  =x[0] W_{8}^{0} + x[1] W_{8}^{2} - x[2] W_{8}^{0}-x[3] W_{8}^{2} + x[4] W_{8}^{0} + x[5] W_{8}^{2}- x[6] W_{8}^{0} - x[7] W_{8}^{2} \qquad$$
$$X[3]  =x[0] W_{8}^{0} + x[1] W_{8}^{3} - x[2] W_{8}^{2}+x[3] W_{8}^{1} - x[4] W_{8}^{0} - x[5] W_{8}^{3}+ x[6] W_{8}^{2} - x[7] W_{8}^{1} \qquad$$
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
The batterfly requires only one complex multiplication $W_N^rX_{m-1}[q]$ and $N\log_2N$ multiplication for computing FFT overall. For example, 8-point FFT requires 24  complex multipliers as can be seen from Fig. (1.2).



<br/>


## Mathematical modelling and SNR explanation




## High Level Synthesis implementation with C/C++





## References

1. source
2. source
3. source