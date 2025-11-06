import os
import matplotlib.pyplot as plt 
from matplotlib.colors import SymLogNorm
import numpy as np 
import pandas as pd
import os

dataPath = "./concatenated_data.csv"
wd = os.getcwd()
outputPath = os.path.join(wd,"PT-analysis-concat-plots")
os.makedirs(outputPath,exist_ok = True)



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


# TODO: Set x and y ticks from the J and D arrays
def plot_image(imageBuff, filename, origin: str = "lower", titleString: str = None, logCmap = False):
    import matplotlib.pyplot as plt 
    from matplotlib.colors import SymLogNorm
    fig, ax = plt.subplots()
     
    if logCmap:
        # We want to fix the colormap for all plots 
        maxVal = np.max(wnd_num)
        minVal = np.min(wnd_num)
        

        vmax = np.max(np.abs([maxVal,minVal]))
        norm = SymLogNorm(linthresh = 0.1*vmax, vmin = minVal, vmax = maxVal)

        mappable = ax.imshow(imageBuff, origin = origin, norm = norm)
    cbar = plt.colorbar(mappable, ax = ax)
    
    ax.set_xlabel("J (ev)")
    ax.set_ylabel("D (ev)")
    fig.savefig(filename,bbox_inches = "tight")
    
    titleString = ""
    if titleString != None:
        ax.set_title(titleString)

    
    ax.clear()
    fig.clf()


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

        sk_num_image = np.zeros((len(JVals_unique),len(DVals_unique)))
        wnd_num_image = np.zeros((len(JVals_unique),len(DVals_unique)))
        
        for j, JVal in enumerate(JVals_unique):
            Jindices, = np.nonzero(JVals_plot == JVal) # Now that we know the indices of all values with a specific Jval we need to figure out where these are for each D value.
            DVals_Jindices = DVals_plot[Jindices] # Now we have all of the values along one J column we need to sort by D so that the elements are in the right order
            sortedD_indices = np.argsort(DVals_Jindices)

            sorted_Jindices = Jindices[sortedD_indices]

            sk_num_Js = sk_num_plot[sorted_Jindices] # Skyrmion Numbers for the specific J value
            sk_num_image[j,:] = sk_num_Js

            wnd_num_Js = wnd_num_plot[sorted_Jindices]
            wnd_num_image[j,:] = wnd_num_Js
                 
         
        # Replace JVals_plot, a length NJxND one dimensional array with a list of indices each element maps to.
        filename_sk_num = os.path.join(wd,f"PT-analysis-concat-plots/sk_num_{t}_{b}.pdf")
        filename_wnd_num = os.path.join(wd,f"PT-analysis-concat-plots/wnd_num_{t}_{b}.pdf")

        skNum_titleString = "Skyrmion Number for $\\mu_b B$ = {Bval:.3f}ev, T = {TVal}K"
        wndNum_titleString = "Winding Number for $\\mu_b B$ = {Bval:.3f}ev, T = {TVal}K"
        plot_image(sk_num_image,filename_sk_num, titleString = skNum_titleString, logCmap = True)
        plot_image(wnd_num_image,filename_wnd_num, titleString = wndNum_titleString, logCmap = True)
        # Use np.where with broadcasting.
        print(f"Shape of dataArray = {np.shape(dataArray)}")

        
        

