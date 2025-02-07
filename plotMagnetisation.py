import numpy as np 
import matplotlib.pyplot as plt 
data = np.loadtxt("Magnetisation4.dat")
x = np.linspace(500,1500,len(data))
plt.plot(x,data)
plt.title(r"$k_B = 0.1$")
plt.xlabel(r"T (K)")
plt.ylabel("Mean spin")
plt.show()
