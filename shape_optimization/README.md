# Running experiments
Required: [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) previously installed

### Build your container image locally
============================

To run the shape optimization experiments for sharp interface idenfication first you will need to build a container with [Firedrake](https://www.firedrakeproject.org/) and [SeismicMesh](https://github.com/krober10nd/SeismicMesh). The shell commands bellow will create a parent folder `project`, clone the source repository inside parent folder and build the singularity container `FS.sif`. The build step will take some time to finish since we are compiling firedrake from the source.    

```shell
mkdir project && cd project
git clone https://github.com/Yuri-Albuquerque/fem_level_set_opt.git
sudo singularity build FS.sif /fem_level_set_opt/shape_optimization/recipe.def 
```

### Runnig experiments
======================

From the parent directory one could run an interactive singularity shell: 

```shell
singularity shell --containall -B ./fem_level_set_opt:/home/<your_username>/data FS.sif
```
Then a shell with `Singularity>` prompt name will appear. Accsses data folder and run the shell script `run.sh` to activate firedrake environment, with
```shell
cd data
source run.sh
```
Finally you could run the experiment with
```shell
mpiexec -n <number_of_shots> python3 /shape_optimization/<forward_or_optimization>/?-<experiment>.py
```
