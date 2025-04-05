import numpy as np
import matplotlib.pyplot as plt

# Number of sample points
N = 10000

# -------------------------------
# Method 1: Cartesian approach (original Fortran: using [0,1])
# -------------------------------
# For a fairer comparison across the entire sphere,
# we sample from [-1, 1] instead.
x1 = np.random.uniform(-1, 1, N)
y1 = np.random.uniform(-1, 1, N)
z1 = np.random.uniform(-1, 1, N)
# Normalize each 3D vector
norm1 = np.sqrt(x1**2 + y1**2 + z1**2)
x1_norm = x1 / norm1
y1_norm = y1 / norm1
z1_norm = z1 / norm1

# -------------------------------
# Method 2: Spherical with uniform theta in [0, π] and phi in [0, 2π]
# -------------------------------
theta = np.random.uniform(0, np.pi, N)
phi = np.random.uniform(0, 2*np.pi, N)
x2 = np.sin(theta) * np.cos(phi)
y2 = np.sin(theta) * np.sin(phi)
z2 = np.cos(theta)

# -------------------------------
# Method 3: Spherical with uniform cos(theta)
# -------------------------------
u = np.random.uniform(0, 1, N)
cos_theta = 2*u - 1        # Uniform in [-1, 1]
theta3 = np.arccos(cos_theta)
phi3 = np.random.uniform(0, 2*np.pi, N)
x3 = np.sin(theta3) * np.cos(phi3)
y3 = np.sin(theta3) * np.sin(phi3)
z3 = cos_theta  # already uniformly distributed

# -------------------------------
# Plotting the results
# -------------------------------
fig = plt.figure(figsize=(18, 12))

# Create 3D scatter plots (top row)
ax1 = fig.add_subplot(2, 3, 1, projection='3d')
ax1.scatter(x1_norm, y1_norm, z1_norm, s=1, color='blue', alpha=0.5)
ax1.set_title('Method 1: Cartesian ([-1,1] sampling) normalized')
ax1.set_xlim([-1, 1]); ax1.set_ylim([-1, 1]); ax1.set_zlim([-1, 1])
ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')

ax2 = fig.add_subplot(2, 3, 2, projection='3d')
ax2.scatter(x2, y2, z2, s=1, color='red', alpha=0.5)
ax2.set_title('Method 2: Spherical (θ ∈ [0,π], φ ∈ [0,2π])')
ax2.set_xlim([-1, 1]); ax2.set_ylim([-1, 1]); ax2.set_zlim([-1, 1])
ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.set_zlabel('Z')

ax3 = fig.add_subplot(2, 3, 3, projection='3d')
ax3.scatter(x3, y3, z3, s=1, color='green', alpha=0.5)
ax3.set_title('Method 3: Spherical with uniform cos(θ)')
ax3.set_xlim([-1, 1]); ax3.set_ylim([-1, 1]); ax3.set_zlim([-1, 1])
ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.set_zlabel('Z')

# Histograms of z-coordinates (bottom row)
bins = 50

ax4 = fig.add_subplot(2, 3, 4)
ax4.hist(z1_norm, bins=bins, density=True, color='blue', alpha=0.7)
ax4.set_title('Z Distribution - Method 1 (Cartesian)')
ax4.set_xlabel('Z'); ax4.set_ylabel('Density')

ax5 = fig.add_subplot(2, 3, 5)
ax5.hist(z2, bins=bins, density=True, color='red', alpha=0.7)
ax5.set_title('Z Distribution - Method 2 (θ uniform)')
ax5.set_xlabel('Z'); ax5.set_ylabel('Density')

ax6 = fig.add_subplot(2, 3, 6)
ax6.hist(z3, bins=bins, density=True, color='green', alpha=0.7)
ax6.set_title('Z Distribution - Method 3 (cos(θ) uniform)')
ax6.set_xlabel('Z'); ax6.set_ylabel('Density')

plt.tight_layout()
plt.show()
