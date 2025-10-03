import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd




def render_single(path,outputPath):
    print(f"Plotting {path}")
    df = pd.read_csv(path)

    x = df['x'].to_numpy()
    y = df['y'].to_numpy()
    z = df['z'].to_numpy()

    Sx = df['Sx'].to_numpy()
    Sy = df['Sy'].to_numpy()
    Sz = df['Sz'].to_numpy()

    x_unique = np.unique(x)
    y_unique = np.unique(y)
    z_unique = np.unique(z)

    maxZ = z_unique[-1]
    minZ = z_unique[0]

    halfwayZ = int((maxZ - minZ) / 2)

    df_slice = df[df['z'] == z_unique[halfwayZ]]


    fig, ax = plt.subplots()

    im = ax.tricontourf(x,y,Sz,levels=100,cmap='jet')
    ax.quiver(x,y,Sx,Sy,Sz,scale=50)
    fig.colorbar(im,label='Sz',ax=ax)
    ax.set_xlabel(r'x $(\AA)$')
    ax.set_ylabel(r'y $(\AA)$')

    if outputPath != None:
        head,tail = os.path.split(path)
        tail = tail.replace(".csv",".pdf")
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

def fortran_E(x, d):
    import math
    if x == 0:
        return f"{0:.{d}f}E+00"
    exp = math.floor(math.log10(abs(x)))
    mantissa = x / (10**exp)
    # Shift mantissa to [0.1,1.0)
    mantissa /= 10
    exp += 1
    return f"{mantissa:.{d}f}E{exp:+03d}"

if __name__ == "__main__":
    dataPath = "./concatenated_data.csv"

    df = pd.read_csv(dataPath,na_values="********")
    df = df.fillna(-1)
    JVals = df['J'].to_numpy()
    DVals = df['D'].to_numpy()
    BVals = df['B'].to_numpy()
    TVals = df['T'].to_numpy()
    winding_numbers = df['winding_number_middle'].to_numpy()
    skyrmion_numbers = df[' skyrmion_number_middle'].to_numpy()
    winding_number_spread = df['winding_number_spread'].to_numpy()

    JVals_unique = np.unique(JVals)
    DVals_unique = np.unique(DVals)
    BVals_unique = np.unique(BVals)


    for i in range(0,len(winding_numbers)):
            if winding_numbers[i] != skyrmion_numbers[i]:
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

                BString = f"{BVals[i]:.{d}E}"
                if BString.startswith("0"):
                    BString = BString[1:]
                elif BString.startswith("-0"):
                    BString = "-" + BString[2:]

                TString = fortran_E(TVals[i],d)

                print(f"{JString}_{DString}_{BString}_{TString}")


    wd = os.getcwd()
    outputPath = os.path.join(wd,"Test-Paper-Visualisations")
    #enumerate_spins(outputPath=outputPath)
