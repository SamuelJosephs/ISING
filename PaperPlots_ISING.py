import os 

print("Hello")


def PlotSpinFile(filePath,visDir=None):
    import matplotlib.pyplot as plt 
    import pandas as pd 
    import numpy as np
    print(f"Plotting {filePath}")
    # Here would be the code to actually plot the file
    df = pd.read_csv(filePath)
    z = np.unique(df['z'].to_numpy())
    minZ = np.min(z)
    maxZ = np.max(z)
    # Select only the middle slice in z
    middleZIndex = int(len(z) / 2)
    df = df[df['z'].to_numpy() == z[int(middleZIndex)]]
    zVal = z[int(middleZIndex)]
    x = df['x'].to_numpy()
    y = df['y'].to_numpy()
    Sx = df['Sx'].to_numpy()
    Sy = df['Sy'].to_numpy()
    Sz = df['Sz'].to_numpy()

    # Now Plot quiver of Sx, Sy, Sz with heatmap of Sz
    #First the heatmap of Sz
    plt.figure(figsize=(10, 8))
    plt.tricontourf(x, y, Sz, levels=100, cmap='jet')
    plt.quiver(x, y, Sx, Sy, Sz, scale=50)
    
    plt.colorbar(label='Sz')
    plt.xlabel('X')
    uniqueX = np.unique(x)
    uniqueY = np.unique(y)
    maxX = np.max(uniqueX)
    minX = np.min(uniqueX)
    maxY = np.max(uniqueY)
    minY = np.min(uniqueY)
    xticks = np.linspace(minX,maxX,10)
    yticks = np.linspace(minY,maxY,10)
    plt.xticks(np.unique(xticks))
    plt.yticks(np.unique(yticks))
    plt.ylabel('Y')
    if (visDir is not None):
        head, tail = os.path.split(filePath)
        HeahHead, HeadTail = os.path.split(head)
        tail = HeadTail + "_" + tail
        visFilePath = os.path.join(visDir, tail.replace(".csv", ".png"))
        print(f"Saving to {visFilePath}")
        plt.savefig(visFilePath,bbox_inches='tight')
    plt.close()




wd = os.getcwd()
with os.scandir('.') as entries:
    for entry in entries:
        if entry.is_dir() and entry.name.startswith("od"):
            dirName = os.path.join(wd, entry.name)
            # Now Need to Scan dirName for frame_Num.csv whith maximum Num 
            maxNum = -1
            maxFilePath = ""
            with os.scandir(dirName) as subEntries:
                for subentry in subEntries:
                    if subentry.is_file() and subentry.name.startswith("frame_") and subentry.name.endswith(".csv"):
                          _, numString = subentry.name.split("_")
                          NumString, _ = numString.split(".")
                          Num = int(NumString)
                          if Num > maxNum:
                              maxNum = Num 
                              maxFilePath = os.path.join(dirName, subentry.name)
            visPath = os.path.join(wd,"paper_visualisations")
            os.makedirs(visPath,exist_ok=True)
            PlotSpinFile(maxFilePath,visDir=visPath)

