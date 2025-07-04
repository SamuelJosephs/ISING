import os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt


directory_path = "./output-dir/"
dfs = []
with os.scandir(directory_path) as it: 
    for entry in it:
        if entry.is_file():
            filename = entry.name
            df = pd.read_csv(directory_path + filename)
            dfs.append(df)
df = pd.concat(dfs,ignore_index = True)
df.columns = df.columns.str.strip()
df.replace("********", 0, inplace=True)
print(df)

df_min_T = df[df["T"] == min(df["T"])]

J_vals_min_T = df_min_T["J"].to_numpy(dtype=np.float64)
D_vals_min_T = df_min_T["D"].to_numpy(dtype=np.float64)
B_vals_min_T = df_min_T["B"].to_numpy(dtype=np.float64)
mz_vals_min_T = df_min_T["sz"].to_numpy(dtype=np.float64)
winding_middle_vals_min_T = df_min_T["winding_number_middle"].to_numpy(dtype=np.float64)

skyrmion_number_min_T = df_min_T["skyrmion_number_middle"].to_numpy()
skyrmion_spread_min_T = df_min_T["winding_number_spread"].to_numpy()


# Filter for values close to 
#skyrmion_number_min_T = np.where(abs(skyrmion_spread_min_T) < 2.1,
#                                        skyrmion_number_min_T,np.zeros(np.shape(skyrmion_number_min_T)))
print(f"Winding numbers: {winding_middle_vals_min_T}")
skyrmion_number_min_T = np.where(abs(skyrmion_spread_min_T) < 2.1,
                                        winding_middle_vals_min_T,np.zeros(np.shape(skyrmion_number_min_T)))
J_values_unique = np.unique(J_vals_min_T)
D_vals_unique = np.unique(D_vals_min_T)

winding_matrix = np.zeros((len(D_vals_unique),len(J_values_unique)))
mz_matrix = np.zeros((len(D_vals_unique),len(J_values_unique)))
for i in range(0,len(skyrmion_number_min_T)):
    D_val = D_vals_min_T[i]
    J_val = J_vals_min_T[i]
    J_Index = np.nonzero(J_values_unique == J_val)[0]
    D_Index = np.nonzero(D_vals_unique == D_val)[0]
    winding_matrix[D_Index,J_Index] = skyrmion_number_min_T[i]
    mz_matrix[D_Index,J_Index] = mz_vals_min_T[i]



fig_sk_hm, axes_sk_hm = plt.subplots()
fig_mz_hm, axes_mz_hm = plt.subplots()

z = axes_sk_hm.pcolormesh(
        J_values_unique,
        D_vals_unique,
        winding_matrix,
        shading="auto",
        cmap="inferno"
        )
l = axes_mz_hm.pcolormesh(J_values_unique,
                          D_vals_unique,
                          mz_matrix,
                          shading="auto",
                          cmap="inferno")
axes_sk_hm.set_xlabel("J")
axes_sk_hm.set_ylabel("D")
axes_mz_hm.set_xlabel("J")
axes_mz_hm.set_ylabel("D")


cbar_sk_hm = fig_sk_hm.colorbar(z,ax=axes_sk_hm)
cbar_sk_hm.set_label("Skyrmion Number")
cbar_mz_hm = fig_mz_hm.colorbar(l,ax=axes_mz_hm)
cbar_mz_hm.set_label("mz")

fig_sk_hm.savefig("./sk_hm.pdf",bbox_inches="tight")
fig_mz_hm.savefig("./mz_hm.pdf",bbox_inches="tight")


plt.show()
