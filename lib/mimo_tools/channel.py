# Define o modelo de canal a ser utilizado
import numpy as np
def scatteringchnmtx(rays, tx_array, rx_array):
    """
    Estima o canal ataves atraves do modelo geometrico narrowband utilizando os dados de ray-tracing, utilizando angulo de elevacao e azimute. Retorna a matriz de canal.
    """
    num_tx = tx_array.size # quantidade de elementos de antena no transmissor
    num_rx = rx_array.size # quantidade de elementos de antena no receptor

    h = np.zeros(num_tx * num_rx).reshape(num_rx, num_tx) # define a matriz de canal com tamanho num_tr por num_rx.
    factor = 1/np.sqrt(num_rx * num_tx)
    #tx_element_spacing = tx_array.element_spacing
    #rx_element_spacing = rx_array.element_spacing
    #tx_wave_length = tx_array.wave_length
    #rx_wave_length = rx_array.wave_length

    if (tx_array.formfactory=="UPA" and rx_array.formfactory=="UPA"):
        print("********************UPA")
        
        for n in range(len(rays)):
            departure_omega_x = 2 * np.pi * tx_array.element_spacing * np.sin(rays[n].departure_theta) * np.cos(rays[n].departure_phi)
            departure_omega_y = 2 * np.pi * tx_array.element_spacing * np.sin(rays[n].departure_theta) * np.sin(rays[n].departure_phi)
    
            arrival_omega_x = 2 * np.pi * rx_array.element_spacing * np.sin(rays[n].arrival_theta) * np.cos(rays[n].arrival_phi)
            arrival_omega_y = 2 * np.pi * rx_array.element_spacing * np.sin(rays[n].arrival_theta) * np.sin(rays[n].arrival_phi)
            
            #factor_departure = (1/np.sqrt(num_tx)) 
            #factor_arrival = (1/np.sqrt(num_rx)) 
            #factor = (np.sqrt(num_rx * num_tx)) * rays[n].received_power
            
            factor_departure = (1/np.sqrt(num_tx)) 
            factor_arrival = (1/np.sqrt(num_rx)) 
            #factor = (1/np.sqrt(num_rx * num_tx)) * rays[n].received_power
            factor = (1/np.sqrt(num_rx * num_tx)) * rays[n].freespace_pathloss
            
            #departure
            departure_vec = np.zeros((1, num_tx)) 
            for m in range(len(tx_array.ura)):
                vecx = np.exp(1j * departure_omega_x * np.arange(len(tx_array.ura[m,:])))
                vecy = np.exp(1j * departure_omega_y * np.arange(len(tx_array.ura[:,m])))
                departure_vec = factor_departure *np.matrix(np.kron(vecy, vecx)) 
    
            #arrival
            arrival_vec = np.zeros((1, num_rx)) 
            for m in range(len(rx_array.ura)):
                vecx = np.exp(1j * arrival_omega_x * np.arange(len(rx_array.ura[m,:])))
                vecy = np.exp(1j * arrival_omega_y * np.arange(len(rx_array.ura[:,m])))
                arrival_vec = factor_arrival * np.matrix(np.kron(vecy, vecx))
          
            h = h + (factor) * (arrival_vec.conj().T * departure_vec) 

    if (tx_array.formfactory=="ULA" and rx_array.formfactory=="ULA"):
        print("********************ULA")
        for n in range(len(rays)):
            departure_omega = [-1j * 2 * np.pi * num_elements * tx_array.element_spacing * np.sin(rays[n].departure_theta) * (1/tx_array.wave_length) for num_elements in range(num_tx)]
            arrival_omega = [-1j * 2 * np.pi * num_elements * rx_array.element_spacing * np.sin(rays[n].arrival_theta) * (1/rx_array.wave_length) for num_elements in range(num_rx)]

            departure_vector = np.exp(departure_omega)
            arrival_vector = np.exp(arrival_omega)

            complex_gain = (-1 * rays[n].freespace_pathloss) * np.exp(-1j * 0)

            prod = np.matmul(arrival_vector, departure_vector.conj().T) 
            h = h + (complex_gain * prod)

    return factor * h

def richscatteringchnmtx(num_tx, num_rx):
    """
    Ergodic channel. Fast, frequence non-selective channel: y_n = H_n x_n + z_n.  
    Narrowband, MIMO channel
    PDF model: Rich Scattering
    Circurly Simmetric Complex Gaussian from: 
         https://www.researchgate.net/post/How_can_I_generate_circularly_symmetric_complex_gaussian_CSCG_noise
    """
    sigma = 1
    h = np.sqrt(sigma/2)*(np.random.randn(num_tx, num_rx) + np.random.randn(num_tx, num_rx)*(1j))

    return h
