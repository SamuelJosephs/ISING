import os
import matplotlib.pyplot as plt 
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
# Problem: J values and D values may not be sorted, need to map the J and D vals to indices in the heatmap_buffer array 

def ValsToIndices(JVals,JVal, DVals, DVal):

    maxJ = np.max(JVals)
    minJ = np.min(JVals)

    maxD = np.max(DVals)
    minD = np.min(DVals)

    numJs = len(JVals) # Number of J bins
    numDs = len(DVals) # Number of D bins 
    
    Jrange = maxJ - minJ 
    Drange = maxD - minD

    # JIndex/(numJs - 1) = (Jval - JMin)/(Jrange)
    
    numJs = len(JVals)
    numDs = len(DVals)
    JIndex = np.floor(((JVal - minJ) / (Jrange)) * (numJs - 1))
    DIndex = np.floor(((DVal - minD) / (Drange)) * (numDs - 1))

    return int(DIndex), int(JIndex)


TVals_unique = np.unique(TVals)
BVals_unique = np.unique(BVals)
DVals_unique = np.unique(DVals)
JVals_unique = np.unique(JVals)

heatmap_buffer = np.zeros((len(DVals_unique),len(JVals_unique)))
windingNumber_buffer = np.zeros((len(DVals_unique),len(JVals_unique)))
windingNumberSpread_buffer = np.zeros((len(DVals_unique),len(JVals_unique)))

for TVal in TVals_unique:
    for BVal in BVals_unique:
                
        for i, J in enumerate(JVals):
            JVal = JVals[i]
            DVal = DVals[i]

            DIndex, JIndex = ValsToIndices(JVals_unique, JVal, DVals_unique, DVal)
            if ((winding_number_spread_array[i] <= 10) and (TVals[i] == TVal) and (BVals[i] == BVal)):
                print(f"i = {i} has skyrmion number {skyrmion_number_middle_array[i]}, DIndex, Jindex = {DIndex}, {JIndex}")
                
                heatmap_buffer[DIndex, JIndex] = skyrmion_number_middle_array[i]
                windingNumber_buffer[DIndex,JIndex] = winding_number_middle_array[i]
            if (TVals[i] == TVal) and (BVals[i] == BVal): # Don't want to filter this by spread
                windingNumberSpread_buffer[DIndex,JIndex] = winding_number_spread_array[i]


        # Now for plotting a graph for each T and B value
        fig, ax = plt.subplots()
        im = ax.imshow(heatmap_buffer)
        fig.colorbar(im,ax=ax)
        fig.savefig(f"SKNUM_T_{TVal}_B_{BVal}.pdf",bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots()
        im = ax.imshow(windingNumber_buffer)
        fig.colorbar(im,ax=ax)
        fig.savefig(f"WINDNUM_T_{TVal}_B_{BVal}.pdf",bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots()
        im = ax.imshow(windingNumberSpread_buffer)
        fig.colorbar(im,ax=ax)
        fig.savefig(f"SKNUM_SPREAD_T_{TVal}_B_{BVal}.pdf",bbox_inches="tight")
        plt.close(fig)



