import numpy as np
import matplotlib.pyplot as plt


def plot_spectrum(fname, skip=1):
    data = np.genfromtxt(fname, skip_header=skip)[:10000, :]
    ct = np.vectorize(complex)(data[...,1], data[...,2])
    yf = np.fft.hfft(ct)
    N = data.shape[0]
    T = (data[1, 0]-data[0, 0])*41.341374575751
    freq = np.fft.fftfreq(2*N-2, d=T)*219474.63*2*np.pi
    #fig, ax = plt.subplots()
    plt.xlim([500, 2500])
    #plt.ylim([0, 60])
    plt.plot(freq, np.flip(yf)*freq*(1.0-np.exp(-1052*np.abs(freq)/219474.63)))

plot_spectrum("examples/ir_spectrum_calculations/heom_mol_ir/corr_func.out", 2)
plot_spectrum("examples/ir_spectrum_calculations/heom_cavity_loss_ir/corr_func.out", 2)
plt.show()
