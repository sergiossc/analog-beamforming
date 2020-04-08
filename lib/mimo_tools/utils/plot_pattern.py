import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from numpy.linalg import svd, norm

def PlotPattern(codebook, training_set):
    #for sample_id, sample in training_set.items():
    for cw_id, cw in codebook.items():
        rate = []
        for sample_id, sample in training_set.items():
            r = np.log2(1 + norm(np.matmul(cw.T, sample)))
            print ('rate: ', r)
            rate.append(r)
        plt.plot(rate, '-', label=str(cw_id))
    plt.ylabel('rate')
    plt.xlabel('channel')
    plt.legend()
    plt.show()

#def PlotPattern(codebook, training_set, sets, num_tx, num_rx):
#    for cw_id, samples_id in sets.items():
#        print('cw_id: ', cw_id)
#        cw = codebook[cw_id]
#        samples = [training_set[sample_id] for sample_id in samples_id]
#        max_rate = 0
#        max_sample = np.zeros((num_tx, num_rx))
#        max_sample = max_sample.T
#        for sample in samples:
#            rate = np.log2(1 + norm(np.matmul(cw.T, sample))) 
#            if rate > max_rate:
#                max_rate = rate
#                max_sample = sample
#        #max_rate, max_sample = max([(np.log2(1 + norm(np.matmul(cw.T, training_set[sample_id]))), training_set[sample_id]) for sample_id in samples_id])
#        print(max_rate)
#        print(max_sample)
#        u, d, vh = svd(max_sample)
#        print('vh.shape: ', vh.shape)
#        print('cw.shape: ', cw.shape)
#        
#        fc = 60*(10**9) # freq de operacao 60 GHz
#        wave_length = c/fc # cumprimento de onda
#        distance = wave_length * 0.5
#        
#        #w1 = 1/2 * (1 + 1j) 
#        #w2 = 1/2 * (1 - 1j)
#        
#        w = cw #np.array([w1, w2])
#        
#        theta = np.arange(-180, 180)
#        
#        
#        gain_db1 = np.zeros(len(theta))
#        #gain_db2 = np.zeros(len(theta))
#        for i in range(len(theta)):
#            t = np.deg2rad(theta[i])
#            a1 = [np.exp(-1j * 2 * np.pi * n * distance * np.sin(t)/wave_length) for n in range(len(w))]
#            #a2 = [np.exp(-1j * 2 * np.pi * n * distance * np.sin(t)/wave_length) for n in range(len(w))]
#            product1 = np.matmul(w.T, a1)
#            #product2 = np.matmul(w.T, a2)
#            #gain_db1[i] = 10 * np.log10((norm(product))/np.abs(np.matmul(w.conj().T, w)))
#            gain_db1[i] = 10 * np.log10((np.abs(product1) ** 2)/np.abs(np.matmul(w.conj().T, w)))
#            #gain_db2[i] = 10 * np.log10((np.abs(product2) ** 2)/np.abs(np.matmul(w.conj().T, w)))
#        
#        plt.plot(theta, gain_db1, '-', label=str(cw_id))
#
#    plt.ylabel('Gain(dB)')
#    plt.xlabel('Angle(deg)')
#    plt.legend()
#    plt.show()
