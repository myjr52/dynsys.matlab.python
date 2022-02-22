## Chapter 1 Solutions

### Exercise 1.1
Run the following lines with "matlab_1_1_free_falling_object.m"

```
figure;
plot(tout,xout(:,1));
ylabel('position [m]');
xlabel('time [s]');

figure;
plot(tout,xout(:,2));
ylabel('velocity [m/s]');
xlabel('time [s]');

figure;
plot(tout,xout(:,3));
ylabel('m(t) [kg]');
xlabel('time [s]');
```

### Exercise 1.2
Run the following lines with "matlab_1_1_free_falling_object.m"
```
figure;
subplot(221);
plot(tout,xout(:,1))
ylabel('position [m]');
xlabel('time [s]');
subplot(222);
plot(tout,xout(:,2))
ylabel('velocity [m/s]');
xlabel('time [s]');
subplot(212);
plot(tout,xout(:,3))
ylabel('m(t) [kg]');
xlabel('time [s]');
```

### Exercise 1.3
Run the following lines with "python_1_1_free_falling_object.py". The following codes are based on "https://towardsdatascience.com/matplotlib-multi-column-row-spanning-layouts-f026eb7c3c27"
```
from matplotlib.gridspec import GridSpec
fig = plt.figure(2)
gs = GridSpec(2,2,figure=fig)

ax1 = fig.add_subplot(gs[0,0])
ax1.plot(tout,xout[0,:])
ax1.set_ylabel('position [m]')
ax1.set_xlabel('time [s]')

ax2 = fig.add_subplot(gs[0,1])
ax2.plot(tout,xout[1,:])
ax2.set_ylabel('velocity [m/s]')
ax2.set_xlabel('time [s]')

ax3 = fig.add_subplot(gs[1,:])
ax3.plot(tout,xout[2,:])
ax3.set_ylabel('m(t) [kg]')
ax3.set_xlabel('time [s]')
```
The above codes are advanced levels in matplotlib. One of the good aspects of python is that the syntax and the commands are intuitive to understand. In the above program, *GridSpec* creates 2x2 tiles in the figure, where *gs[0,0]* indicates the first row and the first column of the tiles and *gs[1,:]* indicates the whole columns of the second row of the tiles.
  
### Exercise 1.4
See the solution chapter of the book

### Exercise 1.5
See the solution chapter of the book

### Exercise 1.6
Check "matlab_1_2_ligand_receptor_interactions.m" and "python_1_2_ligand_receptor_interactions.py"
