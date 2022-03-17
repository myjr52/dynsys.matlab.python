## Chapter 3 Solutions

### Exercise 3.1
See the solution in the book

### Exercise 3.2
```
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
