## Chapter 2 Solutions

### Exercise 2.1
Run "matlab_2_2_dqdt_unit_norm_error.m"

### Exercise 2.2
Run "matlab_2_4_stchastic_process.m". Use a large value for the number of realizations, *N_realize*, to make the pdf closer to the truth.
Each plot commands changes the surface plot as follows:

```
>> shading flat
```

![shading flat!](./figures/ex2_2_01.png "shading flat")

```
>> shading interp
```

![shading interp!](./figures/ex2_2_02.png "shading interp")


```
>> shading faceted
```

![shading faceted!](./figures/ex2_2_03.png "shading faceted")


### Exercise 2.3
Run "python_2_4_stchastic_process.py". 

*rstride* and *cstride* set the number of sampling widths in the row and the column directions of the data, respectively. Values greater than 1 speed up the rendering speed of surface plots. It is useful when plotting big size data.

The following three figures are for (rstride, cstride) equal to (1,1), (10,1) or (1,10), respectively:

![(rstride, cstrid)=(1,1)!](./figures/ex2_3_01.png "(rstride, cstrid)=(1,1)")
![(rstride, cstrid)=(10,1)!](./figures/ex2_3_02.png "(rstride, cstrid)=(10,1)")
![(rstride, cstrid)=(1,10)!](./figures/ex2_3_03.png "(rstride, cstrid)=(1,10)")


### Exercise 2.4
Run "matlab_2_5_gyroscope_simulation.py" or "python_2_5_gyroscope_simulation.py". 

### Exercise 2.5
See the solution chapter of the book
