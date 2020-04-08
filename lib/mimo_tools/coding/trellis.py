import numpy as np
import commpy.channelcoding.convcode as cc
import commpy.channelcoding as conv_encode

memory = np.array([3])


g_matrix = np.array([[ 0o10, 0o13  ]])


trellis = cc.Trellis(memory, g_matrix)

print ("Size of the smallest block of input bits that can be encoded using the convolutional code: ", trellis.k)
print ("Size of the smallest block of output bits generated using the convolutional code: ", trellis.n)
print ("Total number of delay elements needed to implement the convolutional encoder: ", trellis.total_memory)
print ("Number of states in the convolutional code trellis: ", trellis.number_states)
print ("Number of branches from each state in the convolutional code trellis: ", trellis.number_inputs)

trellis.visualize()

#coded = conv_encode([1, 0], trellis)
