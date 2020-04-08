import pandas as pd
import sqlite3
import numpy as np
"""
Le os resultados de ensaios realizados no Wireless Insite a partir do arquivo SQLite e retorna um array de raios (ray object) com as informacoes acessiveis atraves de propriedades de um objeto ray.

"""

class ray():
    def __init__(self, path_id, channel_id, departure_phi, departure_theta, arrival_phi, arrival_theta, received_power, freespace_pathloss):
        self.path_id = path_id
        self.channel_id = channel_id
        self.departure_phi = np.deg2rad(departure_phi)
        self.departure_theta = np.deg2rad(departure_theta)
        self.arrival_phi = np.deg2rad(arrival_phi)
        self.arrival_theta = np.deg2rad(arrival_theta)
        self.received_power = received_power
        self.freespace_pathloss  = freespace_pathloss
	

def read_data(chn_id):
    """
    As informacoes obtidas podem ser selecionadas a partir da edicao dastring sql.
    sql = "SELECT P.path_id, channel_id, received_power, time_of_arrival, departure_phi, departure_theta, arrival_phi, arrival_theta, cmp_e_x_r, cmp_e_x_i, cmp_e_y_r, cmp_e_y_i, cmp_e_z_r, cmp_e_z_i, cmp_h_x_r, cmp_h_x_i, cmp_h_y_r, cmp_h_y_i, cmp_h_z_r, cmp_h_z_i, freespace_path_loss, freespace_path_loss_woa, e_theta_r, e_theta_i, e_phi_r, e_phi_i, cir_phs, cmp_volt_r, cmp_volt_i FROM [path] AS P JOIN path_utd ON P.path_id = path_utd.path_id"
    
    Existem outras tabelas com outras informacoes dos resultados dos ensaios. Pode-se usar o DB Browser for SQLite pra visualizar o banco de dados completo.
    """

    conn = sqlite3.connect("/home/snow/github/land/lib/mimo_tools/data/urban_land.urban_land_x3d.sqlite")
    
    sql = "SELECT P.path_id, channel_id, departure_phi, departure_theta, arrival_phi, arrival_theta, received_power, freespace_path_loss FROM [path] AS P JOIN path_utd ON P.path_id = path_utd.path_id WHERE channel_id="+str(chn_id)+" LIMIT "+str(25)
    
    df = pd.read_sql_query(sql, conn)
    conn.close()

    data_dict = df.to_dict()
    
    num_paths = len(data_dict['path_id'])
    rays = []
    for n in range(num_paths):
        path_id = data_dict['path_id'][n] 
        channel_id = data_dict['channel_id'][n] 
        departure_phi = data_dict['departure_phi'][n] 
        departure_theta = data_dict['departure_theta'][n] 
        arrival_phi = data_dict['arrival_phi'][n]
        arrival_theta = data_dict['arrival_theta'][n] 
        received_power = data_dict['received_power'][n] 
        freespace_pathloss = np.power(10, data_dict['freespace_path_loss'][n]/10)
        r = ray(path_id, channel_id, departure_phi, departure_theta, arrival_phi, arrival_theta, received_power, freespace_pathloss)
        rays.append(r) 
    rays = pwr_normalize(rays)
    return rays

def simurays(chn_id):
    """
    Create random rays in range from 0 to 2PI, for test porposes.
    """
    num_paths = 10
    rays = []
    phi_angle_range = np.arange(0, np.pi, np.pi/32)  # with 5 degrees interval
    theta_angle_range = np.arange(0, 2*np.pi, np.pi/32)  # with 5 degrees interval
    for n in range(num_paths):
        path_id = n 
        channel_id = chn_id 
        departure_theta = np.random.choice(theta_angle_range, 1) 
        arrival_theta = np.random.choice(theta_angle_range, 1) 
        departure_phi = np.random.choice(phi_angle_range, 1) 
        arrival_phi = np.random.choice(phi_angle_range, 1) 
        received_power = np.random.rand() 
        freespace_pathloss = np.random.rand() 
       
        r = ray(path_id, channel_id, departure_phi, departure_theta, arrival_phi, arrival_theta, received_power, freespace_pathloss)
        rays.append(r) 
    rays = pwr_normalize(rays)
    return rays

def pwr_normalize(rays):
    pwr = np.array([ray.freespace_pathloss for ray in rays ])
    magnitude = np.sqrt(np.sum(pwr ** 2))
    for ray in rays:
        ray.freespace_pathloss =  ray.freespace_pathloss / magnitude
    return rays

#for ch in range(1, 11):
#    print ('>>> channel ', ch)
#    rays = read_data(ch)
#    departure_angles = [ray.departure_theta for ray in rays]
#    arrival_angles = [ray.arrival_theta for ray in rays]
#    path_loss  = [ray.freespace_pathloss for ray in rays]
#    print ('# of rays: ', len(departure_angles))
#    print(' Departure angles for each ray: ', np.rad2deg(departure_angles))
#    print ('# of rays: ', len(arrival_angles))
#    print(' Arrival angles for each ray: ', np.rad2deg(arrival_angles))
#    print(' Freespace path loss: ', path_loss)
#    print(' Sum freespace path loss: ', np.sum(path_loss))
#
#
