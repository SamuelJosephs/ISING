import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
df = pd.read_csv("./output.csv")

df_min_T = df[df["T"] == min(df["T"])]
print(df_min_T)

J_vals_min_T = df_min_T["J"].to_numpy()
D_vals_min_T = df_min_T["D"].to_numpy()
skyrmion_number_min_T = df_min_T["skyrmion_number_middle"].to_numpy()
skyrmion_spread_min_T = df_min_T["winding_number_spread"].to_numpy()

# Filter for values close to 
skyrmion_number_min_T = np.where(abs(skyrmion_spread_min_T) < 0.9,
                                        skyrmion_number_min_T,np.zeros(np.shape(skyrmion_number_min_T)))
J_values_unique = np.unique(J_vals_min_T)
D_vals_unique = np.unique(D_vals_min_T)

winding_matrix = np.zeros((len(D_vals_unique),len(J_values_unique)))
for i in range(0,len(skyrmion_number_min_T)):
    D_val = D_vals_min_T[i]
    J_val = J_vals_min_T[i]
    J_Index = np.nonzero(J_values_unique == J_val)[0]
    D_Index = np.nonzero(D_vals_unique == D_val)[0]
    winding_matrix[D_Index,J_Index] = skyrmion_number_min_T[i]

plt.pcolormesh(
        J_values_unique,
        D_vals_unique,
        winding_matrix,
        shading="auto",
        cmap="viridis"
        )
plt.show()
