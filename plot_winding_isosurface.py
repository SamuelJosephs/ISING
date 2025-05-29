import numpy as np 
import pandas as pd 
import plotly.graph_objects as go
data = pd.read_csv("./density_array_3d.csv")
print(data)

iso_value = data["Winding_Density"].min()

iso = go.Isosurface(x = data["x"], y = data["y"], z = data["z"], value = data["Winding_Density"], isomin=iso_value, isomax = iso_value/10,
                    colorbar=dict(title="Winding Density"), colorscale="Viridis")
fig = go.Figure(data=iso)
fig.show()
# Save as a self-contained HTML file
fig.write_html("winding_density_isosurface.html", 
               full_html=True,        # include <html> wrapper
               include_plotlyjs="cdn" # pull Plotly.js from the web (smaller file)
)
