# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
r_bin_edges = np.linspace(0, 5, 32)
theta_bin_edges = np.linspace(0, 2*np.pi, 16)

x_bin_edges = []
y_bin_edges = []

for r in r_bin_edges:
    for theta in theta_bin_edges:
        x_bin_edges.append(r * np.cos(theta))
        y_bin_edges.append(r * np.sin(theta))

# %%
x_bin_edges

# %%
plt.plot(x_bin_edges, y_bin_edges, '.')

# %%
r = np.linspace(0, 5, 100)
theta = np.linspace(0, 2*np.pi, 100)
R, Theta = np.meshgrid(r, theta)
Z = np.exp(-R**2 / (2 * 2.0**2))  # Gaussian with sigma=1

polar_gridpts = np.array([(theta, r) for r in r_bin_edges for theta in theta_bin_edges])

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# fig, ax = plt.subplots(1, 1)
plt.plot(polar_gridpts[:, 0], polar_gridpts[:, 1], '.k')
ax.contourf(Theta, R, Z)
plt.show()



# %%
from scipy import integrate

def integrand(r, theta):
    return np.exp(-r**2/2) * r 

binned_avgs = np.zeros((len(r_bin_edges)-1, len(theta_bin_edges)))

r_bin_centers = np.zeros(len(r_bin_edges)-1)
theta_bin_centers = np.zeros(len(theta_bin_edges))

for i in range(len(r_bin_edges)-1):
    for j in range(len(theta_bin_edges)):
        r_bin_centers[i] = (r_bin_edges[i] + r_bin_edges[i+1])/2.
        theta_bin_centers[j] = (theta_bin_edges[j] + theta_bin_edges[(j+1)%len(theta_bin_edges)])/2.
        result, error = integrate.dblquad(integrand, theta_bin_edges[j], theta_bin_edges[(j+1)%len(theta_bin_edges)], lambda r:r_bin_edges[i], lambda r: r_bin_edges[i+1])

        binned_avgs[i, j] = result

        print(r_bin_edges[i], r_bin_edges[i+1])

theta_bin_centers[-1] = (theta_bin_edges[-1] + theta_bin_edges[0] + 2*np.pi)/2.
binned_avgs[0, -1], error = integrate.dblquad(integrand, theta_bin_edges[-1], theta_bin_edges[0] + 2*np.pi, lambda r:r_bin_edges[0], lambda r: r_bin_edges[1])
binned_avgs[1, -1], error = integrate.dblquad(integrand, theta_bin_edges[-1], theta_bin_edges[0] + 2*np.pi, lambda r:r_bin_edges[1], lambda r: r_bin_edges[2])

# %%
bin_edges = np.array([(theta, r) for r in r_bin_edges for theta in theta_bin_edges])

bin_centers = []

for i in range(len(r_bin_centers)):
    for j in range(len(theta_bin_centers)):
        bin_centers.append([r_bin_centers[i], theta_bin_centers[j]])

bin_centers = np.array(bin_centers)

fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})

ax.plot(bin_edges[:, 0], bin_edges[:, 1], '.k')

ax.scatter(bin_centers[:, 1], bin_centers[:, 0], c=binned_avgs)
plt.show()

# %%
x = np.linspace(-2, 2, 20)
y = np.linspace(-2, 2, 20)

X, Y = np.meshgrid(x, y)

R = np.sqrt(X**2 + Y**2)
THETA = np.arctan2(Y, X)

# %%
R

# %%
fig, axs = plt.subplots(subplot_kw={'projection': 'polar'})
plt.plot(THETA, R, '.')

# %%
# Convert the polar bin centers to cartesian
cart_bin_centers = np.zeros((len(bin_centers), 2))
for i in range(len(bin_centers)):
    cart_bin_centers[i] = [bin_centers[i, 0] * np.cos(bin_centers[i, 1]), bin_centers[i, 0] * np.sin(bin_centers[i, 1])]

# %%
# Associate the binned values in polar coordinates to the Cartesian bin centers

cartesian_vals = np.zeros(len(bin_centers))

for i in range(len(r_bin_centers)):
    for j in range(len(theta_bin_centers)):
        m = i*len(theta_bin_centers) + j
        cartesian_vals[m] = binned_avgs[i, j]

# %%
#Compute the nearest coordinate

nearest_coords = np.zeros((len(x), len(y)))

for i in range(len(x)):
    for j in range(len(y)):
        min_dist = (cart_bin_centers[m, 0] - x[i])**2 + (cart_bin_centers[m, 1] - y[j])**2
        min_m = 0
        for m in range(len(cart_bin_centers)):
            dist = (cart_bin_centers[m, 0] - x[i])**2 + (cart_bin_centers[m, 1] - y[j])**2
            print(dist)
            if dist<min_dist:
                min_dist = dist
                min_m = m

        nearest_coords[i, j] = min_m


# %%
nearest_coords

# %%
cart_interpolated = np.zeros(nearest_coords.shape)

# %%
for i in range(len(x)):
    for j in range(len(x)):
        cart_interpolated[i, j] = cartesian_vals[int(nearest_coords[i, j])]

# %%
cart_interpolated

# %%
x

# %%
plt.pcolormesh(X, Y, cart_interpolated)

# %%



