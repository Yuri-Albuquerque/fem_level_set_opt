# Run with firedrake singularity container


Build your container locally
========================================

To run the shape optimization experiments for sharp interface idenfication first you will need to build a container with firedrake and seismic mesh.  

```shell
sudo singularity build firedrake.sif recipe.def 
```