import os
import matplotlib.pyplot as plt 
from matplotlib.colors import SymLogNorm
import numpy as np 
import pandas as pd

dataPath = "./concatenated_data.csv"

df = pd.read_csv(dataPath,na_values="********")
df = df.fillna(-1)
print(df)
winding_number_middle_array = df["winding_number_middle"].to_numpy()
skyrmion_number_middle_array = df["skyrmion_number_middle"].to_numpy()
winding_number_spread_array = df["winding_number_spread"].to_numpy()

JVals = df["J"].to_numpy()
DVals = df["D"].to_numpy()
BVals = df["B"].to_numpy()
TVals = df["T"].to_numpy()
sk_num = df["skyrmion_number_middle"].to_numpy()
wnd_num = df["winding_number_middle"].to_numpy()
sprd_num = df["winding_number_spread"].to_numpy()
# Problem: J values and D values may not be sorted, need to map the J and D vals to indices in the heatmap_buffer array 

dataBuffer = np.array([BVals, # Each of these arrays has the same length, just very non unique
                       TVals,
                       JVals,
                       DVals,
                       sk_num,
                       wnd_num])
print(f"Shape of dataBuffer = {np.shape(dataBuffer)}")

JVals_unique = np.unique(JVals)
DVals_unique = np.unique(DVals)
BVals_unique = np.unique(BVals)
TVals_unique = np.unique(TVals)

J_indices = [n for n,val in enumerate(JVals_unique)]
D_indices = [n for n,val in enumerate(DVals_unique)]

buffShape = (len(BVals_unique),
             len(TVals_unique),
             len(JVals_unique),
             len(DVals_unique))

heatMapBuffers = np.zeros(buffShape)

for t,TVal in enumerate(TVals_unique):
    for b, BVal in enumerate(BVals_unique):
        mask = (dataBuffer[0] == BVal) & (dataBuffer[1] == TVal)
        print(f"mask = {mask}")
        dataBufferSlice = dataBuffer[:,mask]
        print(f"shape of slice = {np.shape(dataBufferSlice)}")
        
        # dataArray has shape (Jvals, DVals, sk_num, wnd_num)
        dataArray = np.array([dataBufferSlice[i] for i in range(2,np.shape(dataBufferSlice)[0])])

        # Need to extract the arrays to plot 
        JVals_plot = dataArray[0]
        DVals_plot = dataArray[1]
        sk_num_plot = dataArray[2]
        wnd_num_plot = dataArray[3]


        # Replace JVals_plot, a length NJxND one dimensional array with a list of indices each element maps to.

        # Use np.where with broadcasting.
        print(f"Shape of dataArray = {np.shape(dataArray)}")

        fig, ax = plt.subplots()

        ax.imshow()

