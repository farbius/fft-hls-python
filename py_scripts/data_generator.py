import numpy as np
import argparse
import os

Nq      = 15
# RMS quantization level
SNR_q   = 6.02*Nq + 1.76


parser = argparse.ArgumentParser()
parser.add_argument('Npoints' ,default=None, type=int)
parser.add_argument('Nsignals',default=None, type=int)
parser.add_argument('SNR_dB'  ,default=None, type=int)
args = parser.parse_args()

# example
# python data_generator.py 1024 4 40

def main():
    print('<< Generating Input Data for FFT')
    print('<< aleksei.rostov@protonmail.com')
    if args.Npoints is None:
        Npoints = 32
    else:
        Npoints = args.Npoints
        
    if args.Nsignals is None:
        Nsignals = 4
    else:
        Nsignals = args.Nsignals
        
    if args.SNR_dB is None:
        SNR_dB = 0
    else:
        SNR_dB = args.SNR_dB


    Ampl     = 2**15
    x        = np.zeros(Npoints, dtype=complex)
    for k in range(Nsignals):
        nbin = Npoints * np.random.rand() + 1
        x   += 10**((SNR_dB - 2*k)/20) * np.exp(2*1j*np.pi*nbin*np.arange(Npoints)/Npoints) + np.random.randn(Npoints)
    x       *= 10**(-SNR_q/20)
    x_16     = np.round(Ampl*x)
    
    
    np.savetxt('../sim_files/data_re.txt', np.real(x_16),fmt='%d')
    np.savetxt('../sim_files/data_im.txt', np.imag(x_16),fmt='%d')
    
    # saving non-scaled input signal
    np.savetxt('../sim_files/data_nonscl_re.txt', np.real(x),fmt='%f')
    np.savetxt('../sim_files/data_nonscl_im.txt', np.imag(x),fmt='%f')
    
    print('<< Successfully Done')
    
    
if __name__ == "__main__":
    main()
