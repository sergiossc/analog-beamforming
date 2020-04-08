import sys

sys.path.append(r'./lib/mimo_tools/')
sys.path.append(r'./lib/mimo_tools/utils')

from numpy.linalg import svd, matrix_rank, eig, norm
import numpy as np
from scipy.constants import c
from read_data import read_data, simurays
from phased import PartitionedArray
from channel import scatteringchnmtx, richscatteringchnmtx
import matplotlib.pyplot as plt

from scipy.io import savemat # save a python dict in mat file.
import pandas as pd
import uuid

from plot_pattern import PlotPattern

if __name__ == "__main__":
    print("60GHz MmWave MIMO.")
   
    num_tx = 4 # numero de antenas no transmissor
    num_tx_rf = 4 # numero de cadeias de rf no transmissor
   
    num_rx = 4 # numero de antenas no receptor
    num_rx_rf = 4 # numero de cadeias de rf no receptor
   
    num_stream = 1 # nummero de streams a serem transmitidos

    fc = 60*(10**9) # freq de operacao 60 GHz
    wave_length = c/fc # cumprimento de onda

    element_spacing = wave_length/2 # espacamento entre os elementos antena no transmissor /receptor MIMO

    # Uniform Lnear Array
    tx_array = PartitionedArray(num_tx, element_spacing, num_tx_rf, wave_length, "ULA") # configuracao da antena de transmissao
    rx_array = PartitionedArray(num_rx, element_spacing, num_rx_rf, wave_length, "ULA") # configuracao da antena de recepcao
    
    num_iteractions = 1000
    num_samples = 10
    num_clusters = 4 # is the k number of codewords

    # getting samples of channel
    training_set = {}
    for s in range(1, num_samples):
        rays = read_data(s)
        h = scatteringchnmtx(rays, tx_array, rx_array)
        sample_id = uuid.uuid4()
        training_set[sample_id] = h

    # getting the initial random codebook 
    codebook = {}
    np.random.seed(int(sys.argv[1]))
    for k in range(num_clusters):
        cw_id = uuid.uuid4()
        phase_shifters = np.deg2rad(np.random.choice(180, num_rx, replace=False)) # getting num_rx random phase-shifters from 0 to 180 degres then convert to radian 
        codebook[cw_id] = np.exp(1j * phase_shifters) # convert the phase-shifters in complex numbers with unitary modulus

    sum_rate = np.zeros(num_iteractions) # This is que function who I want to get max value
    sets = {}
    for n in range(num_iteractions):
        print('>>> Iteraction #', n)    
        # clustering each of channel sample in a cell 
        for cw_id, cw in codebook.items():
            sets[cw_id] = []

        for sample_id, sample in training_set.items():
            max_rate, cw_id = max([(np.log2(1 + norm(np.matmul(cw.T, sample))), cw_id) for cw_id, cw in codebook.items()])
            sets[cw_id].append(sample_id)

        # Getting sum of metric for measurements porpose 
        sum_max_rate = 0
        for cw_id, samples_id in sets.items():
            cw = codebook[cw_id]
            samples = [training_set[sample_id] for sample_id in samples_id]
            sum_max_rate =+ np.sum([np.log2(1 + norm(np.matmul(cw.T, sample))) for sample in samples])
        sum_rate[n] = sum_max_rate 

        # Getting a new codebook who may improves the rate
        for cw_id, samples_id in sets.items():
            cw = codebook[cw_id]
            current_phase = np.angle(cw)
            samples = [training_set[sample_id] for sample_id in samples_id]
            epison = 0.01
            sum_adjust = 0
            for sample in samples:
                adjust = np.angle((2 * np.matmul(cw.T, np.matmul(sample.T, sample)) / (((norm(np.matmul(cw.T, sample))) * np.log(2)) + np.log(2))) * ( 1j * np.exp([phase * (1j) for phase in current_phase])))
                sum_adjust =+  adjust

            if (len(samples)>0):
                phase_adjust = (1/len(samples)) * sum_adjust
                new_phase = current_phase + (epison * phase_adjust) 
                codebook[cw_id] = np.exp(-1j * new_phase)

    import matplotlib.pyplot as plt
    plt.plot(sum_rate)
    plt.ylabel('sum rate')
    plt.xlabel('# iteractions')
    #plt.legend()
    plt.show()
 
    df = pd.DataFrame(data=codebook)
    df.to_csv('codebook.csv')
    
    #filename = 'codebook.mat'
    #savemat(filename, dict(codebook=codebook))
#PlotPattern(codebook, training_set, sets, num_tx, num_rx) 
