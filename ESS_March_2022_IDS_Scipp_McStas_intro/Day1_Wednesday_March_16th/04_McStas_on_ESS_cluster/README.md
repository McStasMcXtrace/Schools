## Parallelisation, ESS cluster

As you may have seen when performing the earlier simulations, mcgui
makes use of the underlying ```mcrun``` tool, see below screenshot:

[mcgui screenshot](mcgui.png)

Among the output, you find this line
```
mcrun BNL_H8.instr -d /home/jovyan/work/BNL_H8_20220307_094456 -n 1000000  lambda=2.36

```

Our work in this exericise is mainly to get aqcuainted with using the 
mcrun utility that has many options relevant to parallelisation.



### Information resources
* Input parameters for the mcrun tool https://github.com/McStasMcXtrace/McCode/wiki/mcrun
* Confluence-based guide https://confluence.esss.lu.se/pages/viewpage.action?pageId=292160881
* The ```mcstas_mcsub_slurm``` command can be used to write a slurm batchfile

### Tasks
* Using the above information resources, please
  * Re-perform one of the simulations from last section using MPI on
    the quark queue
  * Re-perform one of the simulations from last section using OpenACC on
    the GPU queue
  

