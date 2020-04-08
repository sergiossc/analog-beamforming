import pandas as pd
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from read_data import read_data

def plot_rx():
    conn = sqlite3.connect("/home/snow/github/land/lib/mimo_tools/data/urban_land.urban_land_x3d.sqlite")
    sql = "SELECT * from rx"
    df = pd.read_sql_query(sql, conn)
    conn.close()
    
    data_dict = df.to_dict()
    num_records = len(data_dict['rx_id'])
    rx_list = []
    for n in range(num_records):
        rx = {}
        rx['rx_id'] = data_dict['rx_id'][n]
        rx['x'] = data_dict['x'][n]
        rx['y'] = data_dict['y'][n]
        rx_list.append(rx)

    #print(rx_list)
    return rx_list

def plot_tx():
    conn = sqlite3.connect("/home/snow/github/land/lib/mimo_tools/data/urban_land.urban_land_x3d.sqlite")
    sql = "SELECT * from tx"
    df = pd.read_sql_query(sql, conn)
    conn.close()
    
    data_dict = df.to_dict()
    num_records = len(data_dict['tx_id'])
    tx_list = []
    for n in range(num_records):
        tx = {}
        tx['tx_id'] = data_dict['tx_id'][n]
        tx['x'] = data_dict['x'][n]
        tx['y'] = data_dict['y'][n]
        tx_list.append(tx)

    #print(tx_list)
    return tx_list

def plot_devices():    
    tx_list = plot_tx()
    rx_list = plot_rx()
    
    devices_x = []
    devices_y = []
    devices_label = []
    
    for i in range(len(tx_list)):
        x = tx_list[i]['x']
        y = tx_list[i]['y']
        devices_x.append(x)
        devices_y.append(y)
        devices_label.append('tx '+str(tx_list[i]['tx_id']))
    
    for i in range(len(rx_list)):
        x = rx_list[i]['x']
        y = rx_list[i]['y']
        devices_x.append(x)
        devices_y.append(y)
        devices_label.append('rx '+str(rx_list[i]['rx_id']))
    
    plt.plot(devices_x, devices_y, 'ro')
    for i in range(len(devices_label)):
        plt.text(devices_x[i], devices_y[i], devices_label[i])
    plt.show()

def plot_ray(ch_id):
    plt.axes(projection='polar')
    rays = read_data(ch_id)
    #rads = [ray.departure_theta for ray in rays] # angulo de saida do transmissor
    rads = [ray.arrival_theta for ray in rays] # angulo de chegada no receptor
    for radiam in rads:
        plt.polar(radiam, 1, 'o')
    plt.show()

for ch_id in range(1,11):
    plot_ray(ch_id)
