## Chapter 3 Solutions

### Exercise 3.1
See the solution in the book

### Exercise 3.2
```matlab
row_num = 6; % matrix row size
col_num = 6; % matrix column size
row = [2 4 1 3 4 5 2 4 1 2 3 5 6 2 4 6 4 5];
col = [1 1 2 2 2 2 3 3 4 4 4 4 4 5 5 5 6 6];
a_ij = [2 5 2 1 3 4 1 4 5 3 4 2 4 4 2 1 4 1];
A = sparse(row,col,a_ij,row_num,col_num); %
full(A) % check the sparse matrix in the full format
```

### Exercise 3.3

See [matlab_exercise_3_3.m](../matlab/matlab_exercise_3_3.m)

### Exercise 3.4

See [python_exercise_3_4.py](../python/python_exercise_3_4.py)

### Exercise 3.5

See [matlab_exercise_3_3.m](../matlab/matlab_exercise_3_5.m) or [python_3_4_improve_graph_with_new_samples_along_optimal_path.py](../python/python_3_4_improve_graph_with_new_samples_along_optimal_path.py)

### Exercise 3.6
See the solution in the book

### Exercise 3.7

Check the following variable:
```
a_triangle 
b_triangle
```
at lines 75-78 in [matlab_3_6_uav_target_tracking_worst_scenario.m](../matlab/matlab_3_6_uav_target_tracking_worst_scenario.m)

or

at lines 81-82 in [python_3_6_uav_target_tracking_worst_scenario.py](../python/python_3_6_uav_target_tracking_worst_scenario.py)

### Exercise 3.8

See [matlab_3_6_uav_target_tracking_worst_scenario.m](../matlab/matlab_3_6_uav_target_tracking_worst_scenario.m) or [python_3_6_uav_target_tracking_worst_scenario.py](../python/python_3_6_uav_target_tracking_worst_scenario.py)

### Exercise 3.9

Run the following lines after run [matlab_3_6_uav_target_tracking_worst_scenario.m](../matlab/matlab_3_6_uav_target_tracking_worst_scenario.m) 
```matlab
% check the second derivative if it is convex
dJ2=[eval(diff(dJdux0,ux0)) eval(diff(dJdux0,uy0)); 
     eval(diff(dJduy0,ux0)) eval(diff(dJduy0,uy0))];

ux_all = linspace(-20,20,20);
uy_all = linspace(-20,20,19); 
eig_dJ2 = zeros(length(ux_all),length(uy_all));
for idx=1:length(ux_all)
    idx
    for jdx=1:length(uy_all)
        eig_dJ2(idx,jdx) = min(eig(eval(subs(dJ2,[ux0 uy0],[ux_all(idx) uy_all(jdx)]))));
    end
end
figure; clf; contourf(ux_all,uy_all, eig_dJ2'); colorbar;
```

Run the following lines after run [python_3_6_uav_target_tracking_worst_scenario.py](../python/python_3_6_uav_target_tracking_worst_scenario.py)

```python
# check the second derivative if it is convex
d2Jux02 = diff(dJdux0,ux0).subs(values)
d2Jux02_function = lambdify([ux0,uy0],d2Jux02)

d2Juy02 = diff(dJduy0,uy0).subs(values)
d2Juy02_function = lambdify([ux0,uy0],d2Juy02)

d2Jux0uy0_1 = diff(dJdux0,uy0).subs(values)
d2Jux0uy0_1_function = lambdify([ux0,uy0],d2Jux0uy0_1)

d2Jux0uy0_2 = diff(dJduy0,ux0).subs(values)
d2Jux0uy0_2_function = lambdify([ux0,uy0],d2Jux0uy0_2)

a11 = d2Jux02_function(UX0,UY0) 
a22 = d2Juy02_function(UX0,UY0)
a12_1 = d2Jux0uy0_1_function(UX0,UY0)
a12_2 = d2Jux0uy0_2_function(UX0,UY0)

print('Check if this is zero:',np.sum(np.abs(a12_1-a12_2)))

num_neg_eig = 0
for ix, iy in np.ndindex(a11.shape):
    J_Hessian = np.array([[a11[ix,iy],a12_1[ix,iy]],[a12_1[ix,iy],a22[ix,iy]]])
    [eig_val,eig_vec]=np.linalg.eig(J_Hessian)
    if eig_val.min() < 0:
        num_neg_eig = num_neg_eig + 1
        
if num_neg_eig == 0:
    print('Hessian is positive-definite')
else:
    print('find negative hessian')
```

### Exercise 3.10

See [matlab_3_7_uav_target_tracking_mc_simulations](../matlab/matlab_3_7_uav_target_tracking_mc_simulations.m) or [python_3_7_uav_target_tracking_mc_simulations.py](../python/python_3_7_uav_target_tracking_mc_simulations.py)
