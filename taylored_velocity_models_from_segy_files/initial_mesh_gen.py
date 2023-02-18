# Compute the difference of two SDFs to create more complex geometries.
import os
import pdb
import meshio
import segyio
import h5py as h5
import SeismicMesh
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from mpi4py import MPI
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Original stratigraphy model
fname = './segy_files/sigsbee2b_stratigraphy.segy'
true_background = './true/original/sigsbee2b_true_background.hdf5'
background = np.array(h5.File(true_background, "r")['velocity_model'])

# Function to read segy file copy from 
#(https://github.com/krober10nd/SeismicMesh/blob/03839c15276bc10f4d7cefc01a9c6efa3ba986e9/SeismicMesh/sizing/mesh_size_function.py#L572)
def read_segy(filename):
    """Read a velocity model from a SEG-y file"""
    import segyio

    with segyio.open(filename, ignore_geometry=True) as f:
        f.mmap()
        nz, nx = len(f.samples), len(f.trace)
        vp = np.zeros(shape=(nz, nx))
        for index, trace in enumerate(f.trace):
            vp[:, index] = trace
        if np.amin(vp) < 1000.0:
            warnings.warn(
                "Velocity appear to be in km/s. Maybe pass `units` km-s key pair?"
            )
        return 0.3*np.flipud(vp), nz, nx, 0

# Save salt shape as a reference
segyf, nz, nx, _ = read_segy(fname)
bbox=(-9151.62, 0.0, 0.0, 27439.62)

# removing hard layer from bottom 
# segy_file[:25,:] = segy_file[25:50, :]
segyf[:21,:] = segyf[21, :]
# segy_file = segyf[20:, :]
segy_file = segyf[:, :]

ref_image = np.zeros_like(segy_file)
ref_image = np.where(segy_file < 4000.0, 0, 1)
inv_ref_image = 1.0 - ref_image
im = Image.fromarray(np.uint8(cm.binary(inv_ref_image)*255))
im.save("./true/ref_image.png")
# Then create a new initial guess from the reference
initial_guess_img = 1.0 - \
    np.int64(mpimg.imread('./guess/init_guess_img.png')==0)[:,:,0]
# np_segy_file = np.where(initial_guess_img == 0, segy_file, 4500)
np_segy_file = np.where(initial_guess_img == 0,
                        background[50:, 62:nx+62], 4500)

# Saving the modified background
segyio.tools.from_array('./guess/background_3.segy',
                        np.fliplr(background[50:, 62:nx+62].T))
# Write the velocity model as an hdf5 file
SeismicMesh.sizing.write_velocity_model('./guess/background_3.segy', ofname='./guess/background_3',
        bbox = bbox,
        domain_pad=500,
        pad_style="edge",
        )

# Save the initial-guess numpy array as a segy file
fname_guess = './guess/sigsbee2b_guess3.segy'
segyio.tools.from_array(fname_guess, np.fliplr(np_segy_file.T))

rectangle = SeismicMesh.Rectangle(bbox)
# Generates a sizing function
ef_guess = SeismicMesh.sizing.get_sizing_function_from_segy(fname_guess, bbox,
    # cr_max=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
    dt=0.001, # for thergiven :math:`\Delta t` of 0.001 seconds
    # space_order = 2, # assume quadratic elements :math:`P=2`
    wl = 10,
    freq = 2,
    # hmin = 7.62,
    hmin = 28.25,
    units = "m-s",
    # grade=0.15,
    domain_pad=500,
    pad_style="edge",
)

# View of sizing function (To verify if it is as expected)
SeismicMesh.sizing.plot_sizing_function(ef_guess)

# Write the velocity model as an hdf5 file
SeismicMesh.sizing.write_velocity_model(fname_guess, ofname='./guess/sigsbee2b_guess3',
                                        bbox = bbox,
                                        domain_pad=500,
                                        pad_style="edge",
                                        )

# Extracts points and cells from sizing function to generate the unstructured mesh
points_guess, cells_guess = SeismicMesh.generate_mesh(domain=rectangle, 
                edge_length=ef_guess, 
                )

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creating True models
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_true = './true/sigsbee2b_true3.segy'
segyio.tools.from_array(fname_true, np.fliplr(segy_file.T))
ef_true = SeismicMesh.sizing.get_sizing_function_from_segy(fname_true, bbox,
    # cr_max=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
    dt=0.001, # for the given :math:`\Delta t` of 0.001 seconds
    # space_order = 2, # assume quadratic elements :math:`P=2`
    wl = 10,
    freq = 2,
    # hmin = 7.62,
    hmin = 28.25,
    units = "m-s",
    # grade=0.15, 
    domain_pad=500,
    pad_style="edge",
)

# View of sizing function (To verify if it is as expected)
SeismicMesh.sizing.plot_sizing_function(ef_true)

# Write the velocity model as an hdf5 file
SeismicMesh.sizing.write_velocity_model(fname_true, ofname='./true/sigsbee2b_true3',
                                        bbox = bbox,
                                        domain_pad=500,
                                        pad_style="edge",
                                        )

# Extracts points and cells from sizing function to generate the unstructured mesh
points_true, cells_true = SeismicMesh.generate_mesh(domain=rectangle, 
                edge_length=ef_true, 
                )


# write mesh as '.msh' file
if rank == 0:
    meshio.write(
        "./guess/sigsbee2b_guess3.msh",
        meshio.Mesh(points_guess / 1000, [("triangle", cells_guess)]),
        file_format="gmsh22",
        binary = False,
    )

    meshio.write(
        "./true/sigsbee2b_true3.msh",
        meshio.Mesh(points_true / 1000, [("triangle", cells_true)]),
        file_format="gmsh22",
        binary = False,
    )

print('\n Done')
