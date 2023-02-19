from mpi4py import MPI
import numpy as np
from firedrake import *
import os
import spyro
import gc
import copy
import SeismicMesh
import meshio 
import segyio
import matplotlib.pyplot as plt
gc.disable()

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1"

model = {}
model["parallelism"] = {
    # options: automatic (same number of cores for evey processor), custom, off
    "type": "automatic",
    # only if the user wants a different number of cores for every shot.
    "custom_cores_per_shot": [],
    # input is a list of integers with the length of the number of shots.
}
model["opts"] = {
    "method": "KMV",
    "degree": 2,  # p order
    "dimension": 2,  # dimension
    "quadrature": "KMV",
}
model["mesh"] = {
    "Lz": 1.0,  # depth in km - always positive
    "Lx": 3.0,  # width in km - always positive
    "Ly": 0.0,  # thickness in km - always positive
    "meshfile": "meshes/shape18_guess_vp.msh",
    "initmodel": "velocity_models/shape18_guess_vp.hdf5",
    "truemodel": "velocity_models/shape18_true_vp.hdf5",
    "background": "velocity_models/shape18_back_vp.hdf5",
}
model["PML"] = {
    "status": True,  # true,  # true or false
    #  dirichlet, neumann, non-reflective (outer boundary condition)
    "outer_bc": "non-reflective",
    "damping_type": "polynomial",  # polynomial, hyperbolic, shifted_hyperbolic
    "exponent": 2,
    "cmax": 4.5,  # maximum acoustic wave velocity in pml - km/s
    "R": 0.001,  # theoretical reflection coefficient
    # thickness of the pml in the z-direction (km) - always positive
    "lz": 0.2,
    # thickness of the pml in the x-direction (km) - always positive
    "lx": 0.2,
    # thickness of the pml in the y-direction (km) - always positive
    "ly": 0.0,
}
recvs = spyro.create_transect((-0.01, 0.250), (-0.01, 2.75), 300)
sources = spyro.create_transect((-0.01, 0.250), (-0.01, 2.750), 3)
model["acquisition"] = {
    "source_type": "Ricker",
    "num_sources": len(sources),
    "source_pos": sources,
    "frequency": 2.0,
    "delay": 0.5,
    "amplitude": 1.0,
    "num_receivers": len(recvs),
    "receiver_locations": recvs,
}
model["timeaxis"] = {
    "t0": 0.0,  # initial time for event
    "tf": 1.5,  # final time for event
    "dt": 0.0001,  # timestep size
    "nspool": 100,  # how frequently to output solution to pvds
    "fspool": 20,  # how frequently to save solution to ram
    # "skip": 2,
}
#### end of options ####

VP_1 = 4.5  # inside subdomain to be optimized

def calculate_indicator_from_vp(vp, mesh):
    """Create an indicator function"""

    dgV = FunctionSpace(mesh, "DG", 0)

    # the salt body is 1 everything else is 2

    cond = conditional(vp > (VP_1 - 0.1), 1, 2)

    indicator = Function(dgV, name="indicator").interpolate(cond)

    return indicator


def update_velocity(V, q, vp_background):
    """Update the velocity (material properties)

    based on the indicator function

    """

    sd1 = SubDomainData(q < 1.5)

    # place background field everywhere

    vp_new = Function(V, name="velocity")

    vp_new.assign(vp_background)

    vp_new.interpolate(Constant(VP_1), subset=sd1)

    return vp_new


def create_weighting_function(V, const=100.0, M=5, width=0.1, show=False):
    """Create a weighting function g, which is large near the

    boundary of the domain and a constant smaller value in the interior


    Inputs

    ------

       V: Firedrake.FunctionSpace

    const: the weight function is equal to this constant value, except close to the boundary

    M:   maximum value on the boundary will be M**2

    width:  the decimal fraction of the domain where the weight is > constant

            can be a scalar value or an array of 4 values corresponding to the width on each side.

    show: Visualize the weighting function



    Outputs

    -------

    wei: a Firedrake.Function containing the weights



    """

    if np.isscalar(width) == 1:

        wleft = wright = wtop = wbottom = width

    elif len(width) == 4:

        wleft, wright, wtop, wbottom = width

    else:

        raise ValueError(
            "Width must be either a scalar or contain four values")

    # get coordinates of DoFs

    m = V.ufl_domain()

    W2 = VectorFunctionSpace(m, V.ufl_element())

    coords = interpolate(m.coordinates, W2)

    Z, X = coords.dat.data[:, 0], coords.dat.data[:, 1]

    a0 = np.amin(X) + 0.0

    a1 = np.amax(X) - 0.0

    b0 = np.amin(Z) + 0.0

    b1 = np.amax(Z) - 0.0

    cx = a1 - a0  # x-coordinate of center of rectangle

    cz = b1 - b0  # z-coordinate of center of rectangle

    def h(t, d, _width):

        L = _width * d  # fraction of the domain where the weight is > constant

        return (np.maximum(0.0, M / L * t + M)) ** 2

    w = const * (

        1.0

        + np.maximum(

            h(X - a1, cx, wright) + h(a0 - X, cx, wleft),

            h(b0 - Z, cz, wtop) + h(Z - b1, cz, wbottom),

        )

    )

    if show:

        import matplotlib.pyplot as plt

        plt.scatter(Z, X, 5, c=w)

        plt.colorbar()

        plt.clim([1e-10, 1e-8])

        plt.savefig("weighting_function.png")

    wei = Function(V, w, name="weighting_function")

    File("./results/weighting_function/weighting_function.pvd").write(wei)

    # quit()

    return wei


def calculate_functional(model, mesh, comm, vp, sources, receivers, iter_num):
    """Calculate the l2-norm functional"""

    if comm.ensemble_comm.rank == 0:
        print("Computing the cost functional", flush=True)

    J_local = np.zeros((1))

    J_total = np.zeros((1))

    for sn in range(model["acquisition"]["num_sources"]):

        if spyro.io.is_owner(comm, sn):

            # calculate receivers for each source

            guess, guess_dt, guess_recv = spyro.solvers.Leapfrog_level_set(

                model, mesh, comm, vp, sources, receivers, source_num=sn

            )

            f = "shots/forward_exact_level_set" + str(sn) + ".dat"

            p_exact_recv = spyro.io.load_shots(f)

            # DEBUG

            # viz the signal at receiver # 100

            if comm.comm.rank == 0:

                import matplotlib.pyplot as plt

                plt.plot(p_exact_recv[:-1:1, 100], "k-")

                plt.plot(guess_recv[:, 100], "r-")

                # plt.ylim(np.min(guess_recv[:,100])*(1-0.05), np.max(guess_recv[:,100])*(1+0.05))
                plt.ylim(-0.01, 0.01)

                plt.title("Receiver #100")

                plt.savefig(

                    "./results/comparison_"

                    + str(comm.ensemble_comm.rank)

                    + "_iter_"

                    + str(iter_num)

                    + ".png"

                )

                plt.close()

            # END DEBUG

            residual = spyro.utils.evaluate_misfit(

                model,

                comm,

                guess_recv,

                p_exact_recv,

            )

            J_local[0] += spyro.utils.compute_functional(model, comm, residual)

    if comm.ensemble_comm.size > 1:

        COMM_WORLD.Allreduce(J_local, J_total, op=MPI.SUM)

        J_total[0] /= comm.ensemble_comm.size

    if comm.ensemble_comm.rank == 0:

        print(f"The cost functional is: {J_total[0]}")

        with open("./results/cost_history.txt", "a+") as ch:

            ch.write(f"{J_total[0]:.50f}\n")

    return (

        J_total[0],

        guess,

        guess_dt,

        residual,

    )


def calculate_gradient(

    model, mesh, comm, vp, vp_background, guess, guess_dt, weighting, residual

):
    """Calculate the shape gradient"""

    if comm.ensemble_comm.rank == 0:

        print("Computing the shape derivative...", flush=True)

    VF = VectorFunctionSpace(
        mesh, model["opts"]["method"], model["opts"]["degree"])

    theta = Function(VF, name="gradient")

    for sn in range(model["acquisition"]["num_sources"]):

        if spyro.io.is_owner(comm, sn):

            theta_local = spyro.solvers.Leapfrog_adjoint_level_set(

                model,

                mesh,

                comm,

                vp,

                vp_background,

                guess,

                guess_dt,

                weighting,

                residual,

                source_num=sn,

                output=False,

                piecewise_smooth=True,

            )

    # sum shape gradient if ensemble parallelism here

    if comm.ensemble_comm.size > 1:

        comm.ensemble_comm.Allreduce(

            theta_local.dat.data[:], theta.dat.data[:], op=MPI.SUM

        )

    else:

        theta = theta_local

    # scale factor

    # theta.dat.data[:] *= -1

    return theta


def model_update(mesh, indicator, theta, step, timesteps=10):
    """Solve a transport equation to move the subdomains around based

    on the shape gradient which hopefully minimizes the functional.

    """

    if comm.ensemble_comm.rank == 0:

        print("Updating the shape...", flush=True)

    indicator_new = spyro.solvers.advect(
        mesh,
        indicator,
        step * theta,
        number_of_timesteps=timesteps,
        output=False,
    )
    return indicator_new


def optimization(
    model, mesh, V, comm, vp, vp_background, sources, receivers, max_iter=10
):
    """Optimization with steepest descent using a line search algorithm"""
    factor = 1.0
    beta0 = beta0_init = 1.5*factor
    max_ls = 8
    gamma = gamma2 = 0.8
    timesteps = 10

    model_new = copy.deepcopy(model) # for the re-mesh process
    # sources_new = copy.deepcopy(sources) # for the re-mesh process
    # receivers_new = spyro.Receivers(model, mesh, V, comm).create() # for the re-mesh process

    # the file that contains the shape gradient each iteration
    if comm.ensemble_comm.rank == 0:
        grad_file = File("./results/theta/theta.pvd", comm=comm.comm)

    weighting = create_weighting_function(
        V, width=(0.1, 0.1, 0.1, 0.1), M=10, const=1e-10, show=False)

    ls_iter = 0

    iter_num = 0

    # calculate the new functional for the new model

    J_old, guess, guess_dt, residual = calculate_functional(

        model, mesh, comm, vp, sources, receivers, iter_num
    )

    while iter_num < max_iter:

        # create a new weighting function

        if iter_num > 100:

            factor = 10

            beta0 = beta0_init = 1.5*factor

            weighting = create_weighting_function(
                V, width=(0.1, 0.1, 0.1, 1.6), M=10, const=1e-10)

        if comm.ensemble_comm.rank == 0 and iter_num == 0 and ls_iter == 0:

            print("Commencing the inversion...")

        if comm.ensemble_comm.rank == 0:

            print(f"The step size is: {beta0}", flush=True)

        # compute the shape gradient for the new domain

        if ls_iter == 0:

            theta = calculate_gradient(

                model,

                mesh,

                comm,

                vp,

                vp_background,

                guess,

                guess_dt,

                weighting,

                residual,

            )

        # write the gradient to a vtk file

        if comm.ensemble_comm.rank == 0:

            grad_file.write(theta, name="gradient")
        
        

        # calculate the so-called indicator function by thresholding vp

        indicator = calculate_indicator_from_vp(vp, mesh)

        # update the new shape by solving the transport equation with the indicator field

        indicator_new = model_update(
            mesh, indicator, theta, beta0, timesteps=timesteps)

        # update the velocity according to the new indicator

        vp_new = update_velocity(V, indicator_new, vp_background)

        
        # write ALL velocity updates to a vtk file
        if comm.ensemble_comm.rank == 0:
        
            evolution_of_velocity.write(vp_new, name="velocity")

        # compute the new functional
        J_new, guess_new, guess_dt_new, residual_new = calculate_functional(
            model, mesh, comm, vp_new, sources, receivers, iter_num
        )

        # write the new velocity to a vtk file

        # using a line search to attempt to reduce the functional

        if J_new < J_old:

            if comm.ensemble_comm.rank == 0:

                print(

                    f"Iteration {iter_num}: Functional was {J_old}. Accepting shape update...new functional is: {J_new}",

                    flush=True,
                )
            # remeshed = False
            ###################################
            #### begin of re-mesh process  ####
            ###################################
            if iter_num > 0 and iter_num % 2 == 0:
                print(f'iteration number: {iter_num}')
                _, _, data_new = spyro.io.write_function_to_grid(vp_new, V, 0.01)
                # plt.imshow(data_new)
                # plt.colorbar()
                # plt.show()
                segyio.tools.from_array('./results/evolution/vp_new_'+str(iter_num)+'.segy', np.fliplr(data_new))
                bbox = (-1000*model['mesh']['Lz'],0.0,0.0,1000*model['mesh']['Lx'])
                rectangle = SeismicMesh.Rectangle(bbox)
                #save hdf5 file
                SeismicMesh.sizing.write_velocity_model('./results/evolution/vp_new_'+str(iter_num)+'.segy', 
                                            ofname='./results/evolution/vp_new_'+str(iter_num),
                                            bbox=bbox,
                                            domain_pad=1000*model['PML']['lz'],
                                            pad_style="edge",
                                            comm=comm.ensemble_comm,
                                            )
                
                
                ef_guess = SeismicMesh.sizing.get_sizing_function_from_segy("./results/evolution/vp_new_"+str(iter_num)+".segy", 
                                                                    bbox,
                                                                    # cr_max=0.5, # maximum bounded Courant number :math:`Cr_{max}` in the mesh
                                                                    dt=0.0001,  # for thergiven :math:`\Delta t` of 0.001 seconds
                                                                    # space_order = 2, # assume quadratic elements :math:`P=2`
                                                                    wl=20,
                                                                    freq=2.0,
                                                                    hmin=15.0,
                                                                    # hmin = 28.25,
                                                                    units="m-s",
                                                                    grade=0.1,
                                                                    domain_pad=200.0,
                                                                    pad_style="edge",
                                                                    comm=comm.ensemble_comm,
                                                                    )
                pointsNew, cellsNew = SeismicMesh.generate_mesh(domain=rectangle,
                                                            edge_length=ef_guess,
                                                            comm=comm.ensemble_comm,
                                                            )
                if comm.ensemble_comm.rank == 0:
                    meshio.write("./results/evolution/evolution_of_velocity_"+str(iter_num)+".msh",
                                    meshio.Mesh(pointsNew/1000, 
                                    [('triangle',cellsNew)]),
                                    file_format="gmsh22",
                                    binary=False,
                                )
                comm.ensemble_comm.Barrier()
                model_new["mesh"]["meshfile"] = "./results/evolution/evolution_of_velocity_"+str(iter_num)+".msh" 
                model_new["mesh"]["initmodel"] = "./results/evolution/vp_new_"+str(iter_num)+".hdf5" 
                mesh, V = spyro.io.read_mesh(model_new, comm)
                vp_new = spyro.io.interpolate(model_new, mesh, V, guess=True)
                vp_background = spyro.io.interpolate(model_new, mesh, V, background=True)
                sources = spyro.Sources(model_new, mesh, V, comm).create()
                receivers = spyro.Receivers(model_new, mesh, V, comm).create()
                model = copy.deepcopy(model_new)
                weighting = create_weighting_function(
                                V, width=(0.1, 0.1, 0.1, 1.6), 
                                M=10, 
                                const=1e-10,
                                )
                # # First, grab the mesh.
                # m = V.ufl_domain()

                # # Now make the VectorFunctionSpace corresponding to V.
                # W = VectorFunctionSpace(m, V.ufl_element())

                # # Next, interpolate the coordinates onto the nodes of W.
                # X = interpolate(m.coordinates, W)

                # # Make an output function.
                # vp_new_msh = Function(V)
                # vp_background_new_msh = Function(V)
                
                # # Use the external data function to interpolate the values of f.
                # vp_new_msh.dat.data[:] = vp_new(X.dat.data_ro)
                # vp_background_new_msh.dat.data[:] = vp_background(X.dat.data_ro)
                # remeshed = True

            ###################################
            ####  end of re-mesh process   ####
            ###################################

            iter_num += 1
            # accept new domain

            J_old = J_new

            guess = guess_new

            guess_dt = guess_dt_new

            residual = residual_new

            indicator = indicator_new
            
            vp = vp_new
            # if remeshed:
            #     vp = vp_new_msh
            #     vp_background = vp_background_new_msh
            # else:
            #     vp = vp_new
        
            # update step

            if ls_iter == max_ls:

                beta0 = max(beta0 * gamma2, 0.1 * beta0_init)

            elif ls_iter == 0:

                beta0 = min(beta0 / gamma2, 1.0 * factor)

            else:

                # no change to step

                beta0 = beta0

            ls_iter = 0

        elif ls_iter < max_ls:

            if comm.ensemble_comm.rank == 0:

                print(

                    f"Previous cost functional {J_old}...new cost functional {J_new}",

                    flush=True,

                )

            # advance the line search counter

            ls_iter += 1

            if abs(J_new - J_old) < 1e-16:

                print(

                    f"Line search number {ls_iter}...increasing step size...",

                    flush=True,

                )

                timesteps *= 2

            else:

                print(

                    f"Line search number {ls_iter}...reducing step size...", flush=True

                )

                # reduce step length by gamma

                beta0 *= gamma

            # now solve the transport equation over again

            # but with the reduced step

            # Need to recompute guess_dt since was discarded above

            # compute the new functional (using old velocity field)

            J_old, guess, guess_dt, residual = calculate_functional(

                model, mesh, comm, vp, sources, receivers, iter_num

            )

        else:

            print(

                f"Failed to reduce the functional after {ls_iter} line searches...",

                flush=True,

            )

            break

    return vp


# run the script

comm = spyro.utils.mpi_init(model)


mesh, V = spyro.io.read_mesh(model, comm)


vp = spyro.io.interpolate(model, mesh, V, guess=True)


vp_background = spyro.io.interpolate(model, mesh, V, background=True)


# visualize the updates with this file

if comm.ensemble_comm.rank == 0:

    evolution_of_velocity = File("./results/evolution/evolution_of_velocity.pvd", comm=comm.comm)

    evolution_of_velocity.write(vp, name="velocity")

    # File("vp_background.pvd").write(vp_background)


# Configure the sources and receivers

sources = spyro.Sources(model, mesh, V, comm).create()

receivers = spyro.Receivers(model, mesh, V, comm).create()

# run the optimization based on a line search for max_iter iterations

vp = optimization(

    model, mesh, V, comm, vp, vp_background, sources, receivers, max_iter=150

)
