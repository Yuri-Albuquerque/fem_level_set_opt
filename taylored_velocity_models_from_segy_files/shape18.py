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


def read_segy(filename):
    """Read a velocity model from a SEG-y file"""

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
        return np.flipud(vp), nz, nx, 0


bbox = (-1000.0, 0.0, 0.0, 3000.0)
rectangle = SeismicMesh.Rectangle(bbox)
# Then create a new initial guess from the reference
img = 1.0 - \
    np.int64(mpimg.imread('./true/shape18.png') == 0)[:, :]
true_mat = np.where(img == 0, 2000, 4500)
background = np.ones_like(true_mat)
background *= 2000
nz, nx = true_mat.shape
XX, ZZ = np.meshgrid(np.linspace(0.0, 3.0, nx),
                     np.linspace(0.0, 1.0, nz))
midz, midx = (0.35, 1.5)
ellipse1 = np.sqrt(0.15*np.square((XX-0.75)) + np.square((ZZ - midz))) - 0.2
ellipse2 = np.sqrt(0.15*np.square((XX-2.25)) + np.square((ZZ - midz))) - 0.2
guess_mat = np.where(np.minimum(ellipse1, ellipse2) <=0, 4500, 2000)
path_velocity, path_meshes, path_true, path_guess = (
    '../velocity_models/', '../meshes/', './true/', './guess/')

fname_guess = 'shape18_guess_vp'
fname_true = 'shape18_true_vp'
fname_back = 'shape18_back_vp'
# Saving the true and guessed velocity models as segy files
segyio.tools.from_array(path_true + fname_true +
                        '.segy', np.fliplr(true_mat.T))
segyio.tools.from_array(path_true + fname_back +
                        '.segy', np.fliplr(background.T))
segyio.tools.from_array(path_guess + fname_guess +
                        '.segy', np.fliplr(guess_mat.T))

# Write the velocity models as an hdf5 files
SeismicMesh.sizing.write_velocity_model(path_guess + fname_guess + '.segy', ofname=path_velocity + fname_guess,
                                        bbox=bbox,
                                        domain_pad=200.0,
                                        pad_style="edge",
                                        comm=comm,
                                        )
SeismicMesh.sizing.write_velocity_model(path_true + fname_true + '.segy', ofname=path_velocity + fname_true,
                                        bbox=bbox,
                                        domain_pad=200.0,
                                        pad_style="edge",
                                        comm=comm,
                                        )
SeismicMesh.sizing.write_velocity_model(path_true + fname_back + '.segy', ofname=path_velocity + fname_back,
                                        bbox=bbox,
                                        domain_pad=200.0,
                                        pad_style="edge",
                                        comm=comm,
                                        )
sgy_file, nz, nx, _ = read_segy(path_guess + fname_guess + '.segy')

# Generates a sizing function from segy file
ef_guess = SeismicMesh.sizing.get_sizing_function_from_segy(path_guess + fname_guess + '.segy', bbox,
                                                            # cr_max=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
                                                            dt=0.0001,  # for thergiven :math:`\Delta t` of 0.001 seconds
                                                            space_order = 2, # assume quadratic elements :math:`P=2`
                                                            wl=20,
                                                            freq=2.0,
                                                            hmin=25.0,
                                                            # hmin = 28.25,
                                                            units="m-s",
                                                            grade=0.1,
                                                            domain_pad=200.0,
                                                            pad_style="edge",
                                                            comm=comm,
                                                            )

ef_true = SeismicMesh.sizing.get_sizing_function_from_segy(path_true + fname_true + '.segy', bbox,
                                                           # cr_max=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
                                                           dt=0.0001,  # for thergiven :math:`\Delta t` of 0.001 seconds
                                                           space_order = 2, # assume quadratic elements :math:`P=2`
                                                           wl=20,
                                                           freq=2.0,
                                                           hmin=25.0,
                                                           #    hmin =28.25,
                                                           units="m-s",
                                                           grade=0.1,
                                                           domain_pad=200.0,
                                                           pad_style="edge",
                                                           comm=comm,
                                                           )

# Extracts points and cells from sizing function to generate the unstructured mesh
points_guess, cells_guess = SeismicMesh.generate_mesh(domain=rectangle,
                                                      edge_length=ef_guess,
                                                      comm=comm,
                                                      )
points_true, cells_true = SeismicMesh.generate_mesh(domain=rectangle,
                                                    edge_length=ef_true,
                                                    comm=comm,
                                                    )

# write mesh as '.msh' file
meshio.write(path_meshes + fname_guess + ".msh",
             meshio.Mesh(points_guess / 1000, [("triangle", cells_guess)]),
             file_format="gmsh22",
             binary=False,
             )
meshio.write(path_meshes + fname_true + ".msh",
             meshio.Mesh(points_true / 1000, [("triangle", cells_true)]),
             file_format="gmsh22",
             binary=False,
             )

# View of sizing function (To verify if it is as expected)
SeismicMesh.sizing.plot_sizing_function(ef_true,
                                        comm=comm,
                                        )
SeismicMesh.sizing.plot_sizing_function(ef_guess,
                                        comm=comm,
                                        )

print('\n All done', flush=True)
# if rank == 0:
#     print('\n All done', flush=True)
