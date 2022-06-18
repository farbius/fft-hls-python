import numpy as np
import matplotlib.pyplot as plt


Ampl    = 2**15 - 1



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
    # print(x)
    # print(y_w)
    y_t = x - y_w
    x_t = x + y_w
    # print(x_t/2)
    # print(y_t/2)
    # exit()
    
    return x_t/2, y_t/2



def butter_freq(x, y, w):
    """
    Radix 2 butterfly implementation (decimation-in-time) with rounding to int16 and scaling factor for output
    @param complex x:  FFT complex point (sample)
    @param complex y:  FFT complex point (sample)
    @param complex w:  Complex coefficient
    @return x_t, y_t  complex samples
    """
    y_t = (x - y) * w
    x_t = (x + y) * 1
    return x_t, y_t

## twiddling coefficients
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
            
            
    # print("Re  ")
    # for k in range(Np):
        # print("{:5.0f}".format(np.real(Y_int16[2, k])), end = " ")
    # print(" ")
    # print("Im ")
    # for k in range(Np):
        # print("{:5.0f}".format(np.imag(Y_int16[2, k])), end = " ")
    
    
    # print("""
    
    
    # """)
    
    # print("Re  ")
    # for k in range(Np):
        # print("{:5.0f}".format(np.real(Y_int16[3, k])), end = " ")
    # print(" ")
    # print("Im ")
    # for k in range(Np):
        # print("{:5.0f}".format(np.imag(Y_int16[3, k])), end = " ")
    
    # print("""
    
    
    # """)
    return y



def main():
    print('<< DFT / FFT math modelling')
    print('<< aleksei.rostov@protonmail.com')
    
    hw      = np.loadtxt("../sim_files/dataOUT.txt", dtype=int)
    hw_re   = hw[0::2]
    hw_im   = hw[1::2]
    hw_cmpx = hw_re + 1j*hw_im
    
    x_re    = np.loadtxt("../sim_files/data_re.txt", dtype=int)
    x_im    = np.loadtxt("../sim_files/data_im.txt", dtype=int)
    xcmpx   = x_re + 1j*x_im
    Np = np.size(xcmpx)
    Nstages = int(np.log2(Np))
    Nb      = Nstages
    print("<< Coefficients initialization")
    wk_16 = coef_init(Np)
    print("<< FFT computing")
    print("""
    
    """)
    # py_cmpx  = fft_dit(xcmpx, wk_16, 'reversal') # reversal output order
    py_cmpx  = fft_dit(xcmpx, wk_16, 'normal')     # normal output order
    
    
    # exit()
    print("<< Plotting results")
    plt.figure(num=1, figsize=(10,10))
    
    plt.subplot(211)
    plt.plot(np.abs(py_cmpx), '.-r', label='PYTHON')
    plt.plot(np.abs(hw_cmpx), '.-b', label='HLS')
    plt.legend()
    plt.title("FFT {} points".format(Np))
    plt.xlabel('bin')
    plt.grid()
    

    
    plt.subplot(212)
    plt.plot(np.abs(py_cmpx) - np.abs(hw_cmpx), '.-r')
    plt.title("Error")
    plt.xlabel('bin')
    plt.grid()
    
    plt.tight_layout()
    plt.show()
    print('<< End Modelling')


if __name__ == "__main__":
    main()