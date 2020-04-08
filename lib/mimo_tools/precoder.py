# Define os metodos de precoding/combining
from numpy.linalg import svd, inv, matrix_rank
from waterfilling import waterfilling, channel
import numpy as np

def opt_precoder(h, num_stream):
    """
    Obtem as matrizes de precoding/combining considerando a informacao de canal CSI existente na matriz h. Esta abordagem considera a existencia de CSI tanto no transmissor(CSIT), quanto no receptor(CSIR). Na pratica, no entanto, as caracteristiacas de mmWave ainda nao tornam possivel a obtencao da coerencia de canal entre transmissor e receptor, principalmente em modelos de sistema dinamicos.
    """

    # Realiza a decomposicao em valores singulares do canal para obter as matrizes de precoding/combining.
    u, s, vh = svd(h)
    print ("s matrix from svd(h): ", s)
    #f_opt corresponde a matriz de precoding otima no transmissor.
    beam_direction_opt = vh[:,0:num_stream]
    

    channels = []
    #r = matrix_rank(h)
    for i in range(matrix_rank(h)):
        #ch = channel(s[i], np.random.randn()/100, -100)
        ch = channel(s[i], 1.0e-3, -100)
        channels.append(ch)

    power_opt = waterfilling(channels, 1)

    #w_opt: corresponde a matriz de combining otima no receptor. Utiliza o metodo de combiner mmse.
    #w_opt = (inv((f_opt.conj().T * h.conj() * h.T * f_opt) + (snr * np.eye(num_stream)))) * f_opt.conj().T * h.conj() #combiner mmse
     
    return beam_direction_opt, power_opt # , w_opt

def spatially_sparse_precoding(h, tx_array, rx_array, num_stream, snr):
    pass

def mrt(w, h):
    """
    Maximun Ratio Transmission: uma vez que o transmissor possui conhecimento do canal H e do vetor de combining W, maximizando a energia (o SNR) no receptor, estima-se o precoder F
    """
    pass

def mrc(f, h):
    """
    Maximun Ratio Combiner
    """
    pass
