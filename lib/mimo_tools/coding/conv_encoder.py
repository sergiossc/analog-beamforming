import sys
 
sys.path.append(r'/home/snow/github/land/lib/mimo_tools')
sys.path.append(r'/home/snow/github/land/lib/mimo_tools/utils')


import numpy as np
import commpy.channelcoding.convcode as cc
import commpy.modulation as modulation

from ber import ber_calc
from bit_gen import bit_gen

#N = 10 #number of symbols per the frame
N = 100 #100 #number of symbols per the frame
#message_bits = np.random.randint(0, 2, N) # message
message_bits = bit_gen(1, N) # message

M = 4 # modulation order (QPSK)
k = np.log2(M) #number of bit per modulation symbol
modem = modulation.PSKModem(M) # M-PSK modem initialization 



#Coding
rate = 1/2 # code rate: k/n, k is theoriginal message  number of bits and n is the code word length
L = 7 # length constraint. 
m = np.array([L-1]) # number of delay elements
#generator_matrix = np.array([[0o5, 0o7]]) # generator branches
generator_matrix = np.array([[0o171, 0o133]]) # generator branches
trellis = cc.Trellis(m, generator_matrix)

tb_depth = 5*(m.sum() + 1) # traceback depth


EbNo = 5 # energy per bit to noise power spectral density ratio (in dB)
snrdB = EbNo + 10*np.log10(k*rate) # Signal-to-Noise ratio (in dB)
noiseVar = 10**(-snrdB/10) # noise variance (power)

N_c = 5 # number of trials

BER_soft = np.empty((N_c,))
BER_hard = np.empty((N_c,))
BER_uncoded = np.empty((N_c,))

num_stream = 1

for cntr in range(N_c):
    print("****Trial number: ", cntr) 
    message_bits = bit_gen(num_stream, N) # message
    coded_bits = cc.conv_encode(message_bits, trellis) # encoding
    
    print ("len(message_bits): ", len(message_bits))
    print ("len(coded_bits): ", len(coded_bits))
    
    modulated_uncoded = modem.modulate(message_bits) # modulation (uncoded case)
    modulated = modem.modulate(coded_bits) # modulation
    
    print ("len(modulated_uncoded): ", len(modulated_uncoded))
    print ("len(modulated_coded): ", len(modulated))
#
    Es = np.mean(np.abs(modulated)**2) # symbol energy
    No = Es/((10**(EbNo/10))*np.log2(M)) # noise spectrum density
    print("Es: ", Es)
    print("No: ", No)
#
#
    noisy = modulated + np.sqrt(No/2)*\
        (np.random.randn(modulated.shape[0])+\
         1j*np.random.randn(modulated.shape[0])) # AWGN
#    
    noisy_uncoded = modulated_uncoded + np.sqrt(No/2)*\
        (np.random.randn(modulated_uncoded.shape[0])+\
         1j*np.random.randn(modulated_uncoded.shape[0])) # AWGN (uncoded case)
#
    demodulated_soft = modem.demodulate(noisy, demod_type='soft', noise_var=noiseVar) # demodulation (soft output)
    demodulated_hard = modem.demodulate(noisy, demod_type='hard') # demodulation (hard output)
    demodulated_uncoded = modem.demodulate(noisy_uncoded, demod_type='hard') # demodulation (uncoded case)
#
    decoded_soft = cc.viterbi_decode(demodulated_soft, trellis, tb_depth, decoding_type='unquantized') # decoding (soft decision)
    decoded_hard = cc.viterbi_decode(demodulated_hard, trellis, tb_depth, decoding_type='hard') # decoding (hard decision)
#    print (decoded_soft[:-(L-1)])
#
    NumErr, BER_soft[cntr] = ber_calc(message_bits, decoded_soft[:-(L-1)]) # bit-error ratio (soft decision)
    print ("NumErr: ", NumErr)
    print ("BER_soft: ", BER_soft[cntr])
    NumErr, BER_hard[cntr] = ber_calc(message_bits, decoded_hard[:-(L-1)]) # bit-error ratio (hard decision)
    print ("NumErr: ", NumErr)
    print ("BER_hard: ", BER_hard[cntr])
    NumErr, BER_uncoded[cntr] = ber_calc(message_bits, demodulated_uncoded) # bit-error ratio (uncoded case)
    print ("NumErr: ", NumErr)
    print ("BER_uncoded: ", BER_uncoded[cntr])
#
#mean_BER_soft = np.mean(BER_soft) # averaged bit-error ratio (soft decision)
#mean_BER_hard = np.mean(BER_hard) # averaged bit-error ratio (hard decision)
#mean_BER_uncoded = np.mean(BER_uncoded) # averaged bit-error ratio (uncoded case)
#
#print("Soft decision:\n"+str(mean_BER_soft)+"\n")
#print("Hard decision:\n"+str(mean_BER_hard)+"\n")
#print("Uncoded message:\n"+str(mean_BER_uncoded)+"\n")
