import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functools import partial

def render_single_density(path,outputPath,winding_number,skyrmion_number,outputPath1d,stripe: bool =False):

    print(f"Plotting {path}, outputPath = {outputPath}")

    df = pd.read_csv(path)

    z = df['z'].to_numpy(dtype=np.float64)

    z_unique = np.unique(z)

    maxZ = z_unique[-1]
    minZ = z_unique[0]

    halfwayZ = int(len(z_unique)/2)

    df_slice = df[df['z'] == z_unique[halfwayZ]]

    x = df_slice['x'].to_numpy(dtype=np.float64)
    y = df_slice['y'].to_numpy(dtype=np.float64)
    z = df_slice['z'].to_numpy(dtype=np.float64)


    Winding_Density = df_slice['Winding_Density'].to_numpy(dtype=np.float64)

    #######################################################################


    # Find Location of maximum Skyrmion Number, Take Path through that coordinate and refine as needed
    # We can restric ourselves to plotting a domain where x,y are within a certain radius of the maximum density
    # First: Find Coordinate
    # Second: Create new df_slice filtering by suitable x and y values that are within a specified distance
    maxWinding = np.max(Winding_Density)
    maxWindingIndex = np.nonzero(np.abs(Winding_Density) == np.max(np.abs(Winding_Density)))[0][0]
    
    xVal = x[maxWindingIndex]
    yVal = y[maxWindingIndex]

    df_slice_1d = df_slice[df_slice["y"] == yVal]

    xVals = df_slice_1d["x"].to_numpy(dtype=np.float64)

    windingArray1d = df_slice_1d["Winding_Density"].to_numpy(dtype=np.float64)



    #######################################################################

    fig, ax = plt.subplots()
    

    im = ax.tricontourf(x,y,Winding_Density,levels=100,cmap='jet')
    fig.colorbar(im,label='Winding Density',ax=ax)
    ax.set_xlabel(r'x $(\AA)$')
    ax.set_ylabel(r'y $(\AA)$')
    # Put a red line through where we are taking the 1d paths 
    if stripe:
        yVals = np.zeros(len(xVals))
        yVals[0:] = yVal
        ax.plot(xVals,yVals,color="red")

    head,tail = os.path.split(path)
    tail = tail.replace(".csv","")
    tail = tail + f"_{int(winding_number)}_{int(skyrmion_number)}" + ".pdf"
    #tail = tail.replace(".csv",".pdf")
    savePath = os.path.join(outputPath,tail)

    fig.savefig(savePath,bbox_inches="tight")
    
    
    fig.clf()
    plt.close(fig)
    ax.clear()

    # Now plot 1d cross section for every y along x for the middle z 
    # Use df_slice
 
    
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    z_unique = np.unique(z)

    halfWayX = x_unique[int(len(x_unique)/2)]
    halfWayY = y_unique[int(len(y_unique)/2)]
    halfWayZ = z_unique[int(len(z_unique)/2)]

    # Now we have need to plot windingArray1d against xVals
    # First Figure out where we are going to store it

    head,tail = os.path.split(path)
    tail = tail.replace(".csv","")
    tail = tail + f"_{int(winding_number)}_{int(skyrmion_number)}" + ".pdf"   
    
    savePath = os.path.join(outputPath1d,tail)


    fig, ax = plt.subplots()
    
    ax.plot(xVals,windingArray1d,label = r"Winding Density $q(\vec{r})$")
    fig.savefig(savePath,bbox_inches = "tight")
    fig.clf()
    ax.clear()
    plt.close(fig)



       


def render_single(path,outputPath,winding_number,skyrmion_number):
    print(f"Plotting {path}")
    df = pd.read_csv(path)

    z = df['z'].to_numpy(dtype=np.float64)

    z_unique = np.unique(z)

    maxZ = z_unique[-1]
    minZ = z_unique[0]

    halfwayZ = int((maxZ - minZ) / 2)
    df_slice = df[df['z'] == z_unique[halfwayZ]]

    x = df_slice['x'].to_numpy(dtype=np.float64)
    y = df_slice['y'].to_numpy(dtype=np.float64)
    z = df_slice['z'].to_numpy(dtype=np.float64)

    Sx = df_slice['Sx'].to_numpy(dtype=np.float64)
    Sy = df_slice['Sy'].to_numpy(dtype=np.float64)
    Sz = df_slice['Sz'].to_numpy(dtype=np.float64)


    fig, ax = plt.subplots()

    im = ax.tricontourf(x,y,Sz,levels=100,cmap='jet')
    ax.quiver(x,y,Sx,Sy,Sz,scale=50)
    fig.colorbar(im,label='Sz',ax=ax)
    ax.set_xlabel(r'x $(\AA)$')
    ax.set_ylabel(r'y $(\AA)$')

    head,tail = os.path.split(path)
    tail = tail.replace(".csv","")
    tail = tail + f"_{int(winding_number)}_{int(skyrmion_number)}" + ".pdf"
    #tail = tail.replace(".csv",".pdf")
    savePath = os.path.join(outputPath,tail)

    fig.savefig(savePath,bbox_inches="tight")

    fig.clf()
    plt.close(fig)



def enumerate_spins(outputPath=None):
    from itertools import repeat
    from concurrent.futures import ThreadPoolExecutor
    wd = os.getcwd()

    spinDir = os.path.join(wd,"output_dir-spins")

    pathList = []
    with os.scandir(spinDir) as fileIter:
        for f in fileIter:
            if f.name.endswith(".csv"):
                pathList.append(f.path)


    if outputPath != None:
        os.makedirs(outputPath,exist_ok=True)
    with ThreadPoolExecutor(max_workers=10) as tp:
        tp.map(render_single,pathList,repeat(outputPath))

def fortran_E(x, d: int):
    import math
    from decimal import Decimal, ROUND_HALF_DOWN

    if x == 0:
        return f"{0:.{d}f}E+00"
    exp = math.floor(math.log10(abs(x)))
    mantissa = x / (10**exp)
    # Shift mantissa to [0.1,1.0)
    mantissa /= 10
    #mantissa = round(mantissa,d+1)
    zeros = "0"*int(d)
    precision = "0." + zeros + "1"
    mantissa = Decimal(mantissa).quantize(Decimal(precision),rounding=ROUND_HALF_DOWN)
    exp += 1
    signString = ""
    if x < 0.0:
        signString = "-"
    if exp == 0:
        return f"{signString}{mantissa:.{d}f}"
    else:
        return f"{signString}{mantissa:.{d}f}E{exp:+d}"

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    from itertools import repeat
    dataPath = "./concatenated_data.csv"

    df = pd.read_csv(dataPath,na_values="********")
    #df = df.fillna(-1)
    JVals = df['J'].to_numpy(dtype=np.float64)
    DVals = df['D'].to_numpy(dtype=np.float64)
    BVals = df['B'].to_numpy(dtype=np.float64) 
    TVals = df['T'].to_numpy(dtype=np.float64)
    winding_numbers = df['winding_number_middle'].to_numpy(dtype=np.float64)
    skyrmion_numbers = df['skyrmion_number_middle'].to_numpy(dtype=np.float64)
    winding_number_spread = df['winding_number_spread'].to_numpy(dtype=np.float64)

    JVals_unique = np.unique(JVals)
    DVals_unique = np.unique(DVals)
    BVals_unique = np.unique(BVals)
    TVals_unique = np.unique(TVals)
 
    wd = os.getcwd()
    outputPath = os.path.join(wd,"Test-Paper-Visualisations")   
    outputPath_density = os.path.join(wd,"Test-Paper-Visualisations-Density")
    outputPath_density1d = os.path.join(wd,"Test-Paper-Visualisations-1d-Density")

    spin_dirname = os.path.join(wd,"output-dir_spins")
    density_dirname = os.path.join(wd,"output-dir_density")

    os.makedirs(outputPath,exist_ok=True)
    os.makedirs(outputPath_density,exist_ok=True)
    os.makedirs(outputPath_density1d,exist_ok=True)

    paths = []
    paths_density = []
    windNums = []
    skyrmNums = []
    for i in range(0,len(winding_numbers)):
            if np.isnan(skyrmion_numbers[i]) or np.isnan(winding_numbers[i]):
                print(f"NaN Found: Skyrmion Number: {skyrmion_numbers[i]}, winding_number: {winding_numbers[i]}")
                continue 
           # elif (abs(int(skyrmion_numbers[i])) == abs(int(winding_numbers[i]))):
           #     print(f"Matching SK and WND numbers: {skyrmion_numbers[i]} , {winding_numbers[i]}")

           # elif abs(int(winding_numbers[i])) != abs(int(skyrmion_numbers[i])):
           #     print(f"Matching SK and WND numbers: {skyrmion_numbers[i]} , {winding_numbers[i]}")    

            elif skyrmion_numbers[i] != winding_numbers[i]:
           # elif abs(skyrmion_numbers[i]) > 0.9*np.max(np.abs(skyrmion_numbers)):  
                print(f"Non Matching SK and WND numbers: {int(skyrmion_numbers[i])} , {int(winding_numbers[i])}")
                # Need to construct the file name 
                d = 4
                JString = f"{JVals[i]:.{d}f}"
                if JString.startswith("0"):
                    JString = JString[1:]
                elif JString.startswith("-0"):
                    JString = "-" + JString[2:]
                DString = f"{DVals[i]:.{d}f}"
                if DString.startswith("0"):
                    DString = DString[1:]
                elif DString.startswith("-0"):
                    DString = "-" + DString[2:]

                BString = f"{BVals[i]:.{d}f}"
                if BString.startswith("0"):
                    BString = BString[1:]
                elif BString.startswith("-0"):
                    BString = "-" + BString[2:]

                TString = fortran_E(TVals[i],d)
                
                TString_Index = np.where(TVals_unique == TVals[i])[0][0] + 1
                nameString = f"spins_{JString}_{DString}_{BString}_{TString}.csv"
                nameString_density = f"density_{JString}_{DString}_{BString}_{TString}.csv"
                filename = os.path.join(spin_dirname,nameString)
               
                filename_density = os.path.join(density_dirname,nameString_density)
                if os.path.isfile(filename):
                    print(f"Found File {filename}")
                    paths.append(filename)
                    windNums.append(int(winding_numbers[i]))
                    skyrmNums.append(int(skyrmion_numbers[i]))
                else:
                    print(f"Failed to find file {filename}")
                
                if os.path.isfile(filename_density):
                    print(f"Found File {filename_density}")
                    paths_density.append(filename_density)

    print(f"paths_density = {paths_density}") 



        
        

    
    render_single_density_stripe = partial(render_single_density,stripe=True)
    with ProcessPoolExecutor(max_workers=10) as pp:
        it1 = pp.map(render_single_density_stripe,paths_density,repeat(outputPath_density),windNums,skyrmNums,repeat(outputPath_density1d))
        it2 = pp.map(render_single,paths,repeat(outputPath),windNums,skyrmNums)

        for process in it1:
            if process != None:
                print("Got: ", process)

        for process in it2:
            if process != None:
                print("Got: ", process)

