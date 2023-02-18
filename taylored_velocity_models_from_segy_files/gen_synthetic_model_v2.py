import pdb
import os
import meshio
import segyio
import h5py as h5
import SeismicMesh
import numpy as np
# from PIL import Image, ImageDraw, ImageFont
from mpi4py import MPI
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
os.environ["OMP_NUM_THREADS"] = "1"

# Original stratigraphy model
fname = './segy_files/sigsbee2b_stratigraphy.segy'
true_background = './true/original/sigsbee2b_true_background.hdf5'
background = np.array(h5.File(true_background, "r")['velocity_model'])

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


bbox=(-9151.62, 0.0, 0.0, 27439.62)
rectangle = SeismicMesh.Rectangle(bbox)
# Then create a new initial guess from the reference
img = 1.0 - \
    np.int64(mpimg.imread('./true/true_sigsbee_v2.png') == 0)[:, :, 0]
nz, nx  = img.shape
XX, ZZ = np.meshgrid(np.linspace(0.0, 27.43962, nx),
                     np.linspace(-9.15162,0.0, nz))
true_mat = np.where(img == 0, background[50:, 62:nx+62], 4500)
midz, midx = (-9.915162 - 0) / 4.0, (27.43962) / 2.0
ellipse = np.sqrt(0.05*np.square((XX-midx)) + np.square((ZZ - midz))) - 1.5
guess_mat = np.where(ellipse <= 0, 4500, background[50:, 62:nx+62])
path_velocity, path_meshes, path_true, path_guess = (
    '../velocity_models/', '../meshes/', './true/', './guess/')

pdb.set_trace()
fname_guess = 'small_sigsbee_guess'
fname_true = 'small_sigsbee_true'
fname_back = 'small_sigsbee_back'
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
                                                            # space_order = 2, # assume quadratic elements :math:`P=2`
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
                                                           # space_order = 2, # assume quadratic elements :math:`P=2`
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
