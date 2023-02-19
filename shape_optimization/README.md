# Run with firedrake singularity container
Required: [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) previously installed

Build your container locally
============================

To run the shape optimization experiments for sharp interface idenfication first you will need to build a container with firedrake and seismic mesh.  

```shell
sudo singularity build firedrake.sif recipe.def 
```

Runnig experiments
======================
```shell
singularity run --containall -B /fem_level_set_opt:/data firedrake.sif <your_experiment>.py
```