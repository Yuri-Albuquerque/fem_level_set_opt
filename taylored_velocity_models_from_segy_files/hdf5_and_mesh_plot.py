import meshio
import h5py as h5
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont


# h5file  = '../../original_meshes/sigsbee2b_guess.hdf5'
# mshfile = '../../original_meshes/sigsbee2b_guess.msh'
# h5file = './true/sigsbee2b_true3.hdf5'
# h5file = './true/immersed_disk_back_vp.hdf5'
h5file = '../results/evolution/vp_new.hdf5'
mshfile = '../results/evolution/evolution_of_velocity.msh'
# mshfile = './guess/sigsbee2b_guess3.msh'
# h5file = './true/original/sigsbee2b_true_background.hdf5'
# mshfile = './../../original_meshes/sigsbee2b_true.msh'

#plot hdf5 file
image = h5.File(h5file, "r")
# print(np.array(image)) #to find the key 'velocity_model'
dataset = image['velocity_model']
img = np.array(dataset)

# Plot hdf5 file
plt.figure(figsize=(10,6))
plt.title(f'This {img.shape[0]}x{img.shape[1]} matrix represents the velocity model.')
plt.imshow(img[::-1, :])
plt.colorbar()
plt.grid(color='r', linestyle='-', linewidth=0.25)
plt.show()
plt.close()


# blurring hdf5 file (not a necessary path)
# kernel = np.outer(signal.gaussian(70, 8), signal.gaussian(70, 8))
# img2 = signal.convolve2d(img, kernel, boundary='symm', mode='same')
# plt.figure(figsize=(10,6))
# plt.title(f'This {img2.shape[0]}x{img.shape[1]} matrix represents the velocity model.')
# plt.imshow(img2[::-1, :], cmap='rainbow')
# plt.colorbar()
# plt.show()
# plt.close()

# Plot mesh file 
meshfile = meshio.read(mshfile)
print(f"Mesh Object {meshfile}")
ncells = len(meshfile.cells_dict['triangle'])
print(f"Cells {ncells}")
print(f"Points {len(meshfile.points)}")

triangles = np.array(meshfile.cells_dict['triangle'])
triangles = triangles.reshape(ncells, 3)
x = meshfile.points

triangles[3:8]
plt.figure(figsize=(10,6))
plt.title(f'This is the mesh generated from the velocity model.')
plt.triplot(x[:,1], x[:,0], triangles)
plt.show()
plt.close()
