import numpy as np
import matplotlib.pyplot as plt


Npoints = 1024
nbin    = Npoints / 4 + .9
SNR_dB  = 90

Nq      = 15
SNR_q   = 6.02*Nq + 1.76
print('RMS quantization is {} dB'.format(SNR_q))
SNR_fft = 10*np.log10(Npoints)
SNR_full= SNR_q + SNR_fft

Ampl    = 2**15

def revBits(n, no_of_bits):
    result = 0
    for i in range(no_of_bits):
        result <<= 1
        result |= n & 1
        n >>= 1
    return result

def butter_time_int16(x, y, w):
    """
    Radix 2 butterfly implementation (decimation-in-time) with rounding to int16 and scaling factor for output
    @param complex int16 x:  FFT complex point (sample)
    @param complex int16 y:  FFT complex point (sample)
    @param complex int16 w:  Complex coefficient
    @return x_t, y_t complex  int16 samples
    """
    y_w = np.round((y * w)/Ampl) # rounding back to 16 bits after multiplication
    y_t = x - y_w
    x_t = x + y_w
    return x_t/2, y_t/2

def coef_init(Npoints):
    """
    Twiddling coefficients generation 
    @param int Npoints:  length of FFT
    @return wk_16 complex int16 coefficient array
    """
    wk_16 = np.zeros(Npoints, dtype='complex')
    for k in range(Npoints):
        wk_16[k] = np.round(Ampl * np.exp(-1j*2*np.pi*k / Npoints))
    return wk_16


def fft_dit(x, w, d_order='normal'):
    """
    FFT decimation-in-time implementation 
    @param  complex int16 x: input data
    @param  complex int16 w: twiddling coefficients
    @return complex int16 y: output data
    """
    Np = np.size(x)
    Ns = int(np.log2(Np))
    Y_int16 = np.zeros((Ns + 1, Np), dtype=complex)
    Y_int16[0, :] = x[:]
    for casc in range(Ns):
        d = 0
        for k in range(Np // 2):
            idx_w = int(np.mod(k, 2**casc))*2**(Ns - 1 - casc)
            if np.mod(k, 2**casc) == 0:
                d = 2*k
            idx_1 = revBits(d, Ns)
            idx_2 = revBits(d + 2**casc, Ns)
            d = d + 1
            Y_int16[casc + 1, idx_1], Y_int16[casc + 1, idx_2] = butter_time_int16(Y_int16[casc, idx_1], Y_int16[casc, idx_2], w[idx_w])
            # if(casc == 2):
                # print(w[idx_w])
                # print(idx_w)
    y = np.zeros(Np, dtype=complex)
    if(d_order=='normal'):
        for k in range(Np):
            y[k] = np.round(Y_int16[Ns, revBits(k, Ns)])
    else:
        for k in range(Np):
            y[k] = np.round(Y_int16[Ns, k])
            
    return y

s        = 10**(SNR_dB/20) *  np.exp(2*1j*np.pi*nbin*np.arange(Npoints)/Npoints) + np.random.randn(Npoints)
s       *= 10**(-SNR_q/20)
# n        =   # np.linalg.norm(n)
xn       = s
x_16     = np.round(Ampl*xn)

xf = np.fft.fft(s)/Npoints

wk_16 = coef_init(Npoints)
py_cmpx  = fft_dit(x_16, wk_16, 'normal')/Ampl

plt.figure(num=3,figsize=(10,10))
plt.subplot(211)
plt.plot(np.real(x_16), '.-b')
plt.title('Real part of the scaled input sugnal')
plt.grid()

plt.subplot(212)
plt.plot(np.imag(x_16), '.-r')
plt.title('Imag part of the scaled input sugnal')
plt.grid()


plt.figure(num=2,figsize=(10,10))
plt.subplot(211)
plt.plot(np.real(s), '.-b')
plt.title('Real part of the input sugnal')
plt.grid()

plt.subplot(212)
plt.plot(np.imag(s), '.-r')
plt.title('Image part of the input sugnal')
plt.grid()

plt.figure(num=1,figsize=(10,10))
plt.subplot(211)
plt.plot(20*np.log10(np.abs(xf)), '.-b')
plt.plot(-SNR_full * np.ones(Npoints), '.-r', label='FFT SNR')
plt.plot(-SNR_q * np.ones(Npoints), '.-g', label='RMS quantization level')
plt.title('SW FFT SNR')
plt.ylabel('power, dB')
plt.legend()

plt.grid()

plt.subplot(212)
plt.plot(20*np.log10(np.abs(py_cmpx)), '.-b')
plt.plot(-SNR_full * np.ones(Npoints), '.-r', label='FFT SNR')
plt.plot(-SNR_q * np.ones(Npoints), '.-g', label='RMS quantization level')
plt.title('HW FFT SNR')
plt.ylabel('power, dB')
plt.legend()
plt.grid()

plt.show()