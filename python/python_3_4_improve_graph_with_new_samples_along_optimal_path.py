#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 07:24:28 2021

@author: jongrae
"""

import numpy as np
import matplotlib.pyplot as plt

# number of samples
num_sample = 500

# map size
map_width = 10
map_height = 5

# x,y coordinates of start and destination of the path to be calculated
xy_start = np.array([0,0])
xy_dest = np.array([9,4])

# spread num_sample random points over the map area
xy_points = np.random.rand(num_sample,2)
xy_points[::,0] = xy_points[::,0]*map_width
xy_points[::,1] = xy_points[::,1]*map_height

# stacking them all together with start and destination
xy_points = np.vstack((xy_start,xy_dest,xy_points))
start_node = 0
end_node = 1

# circular obstacle at [3,3], radius 1.5 & define the boundary
obs_xy = [3,3]
obs_rad = 1.5
th = np.arange(0,2*np.pi+0.01,0.01)
x_obs_0 = obs_rad*np.cos(th)+obs_xy[0]
y_obs_0 = obs_rad*np.sin(th)+obs_xy[1]
xy_obs_0 = np.vstack((x_obs_0,y_obs_0)).T

# non-convex obstacle boundary
x_obs_1 = np.array([6,8,8,5,5,7,7,6,6])
y_obs_1 = np.array([1,1,4,4,3,3,2,2,1])
xy_obs_1 = np.vstack((x_obs_1,y_obs_1)).T

# define obstacle using Path in matplotlib.path
from matplotlib.path import Path
Obs_0 = Path(xy_obs_0)
Obs_1 = Path(xy_obs_1)

# found points are not inside the circular obstacle
mask_0 = ~Obs_0.contains_points(xy_points)
xy_points = xy_points[mask_0,::]

mask_1 = ~Obs_1.contains_points(xy_points)
xy_points = xy_points[mask_1,::]

# construct graph using delaunay
from scipy.spatial import Delaunay
tri = Delaunay(xy_points)

# found triangle definition index
temp_idx=tri.simplices[::,0]
temp_jdx=tri.simplices[::,1]
temp_kdx=tri.simplices[::,2]

# remove longer paths, which are likely passing through the obstacle
dist_ij = np.sqrt(np.sum((xy_points[temp_idx,::]-xy_points[temp_jdx,::])**2,1))
dist_jk = np.sqrt(np.sum((xy_points[temp_jdx,::]-xy_points[temp_kdx,::])**2,1))
dist_ki = np.sqrt(np.sum((xy_points[temp_kdx,::]-xy_points[temp_idx,::])**2,1))
dd_all = np.hstack((dist_ij,dist_jk,dist_ki))
cut_dist = np.mean(dd_all)+1*np.std(dd_all)

# distance thereshold for removing longer paths
cut_mask_ij = dist_ij<cut_dist
cut_mask_jk = dist_jk<cut_dist
cut_mask_ki = dist_ki<cut_dist
temp_xy_ij = np.vstack((temp_idx[cut_mask_ij],temp_jdx[cut_mask_ij]))
temp_xy_jk = np.vstack((temp_jdx[cut_mask_jk],temp_kdx[cut_mask_jk]))
temp_xy_ki = np.vstack((temp_kdx[cut_mask_ki],temp_idx[cut_mask_ki]))

# corresponding distance to the paths
dist_ij = dist_ij[cut_mask_ij]
dist_jk = dist_jk[cut_mask_jk]
dist_ki = dist_ki[cut_mask_ki]

# change format into row, column and the distance
xy_index = np.hstack((temp_xy_ij,temp_xy_jk,temp_xy_ki)).T
row_org = xy_index[::,0]
col_org = xy_index[::,1]
row = np.hstack((row_org,col_org))
col = np.hstack((col_org,row_org))
dist = np.hstack((dist_ij,dist_jk,dist_ki))
dist = np.hstack((dist,dist))
num_node = xy_points.shape[0]

# construct the distance matrix
from scipy.sparse import csr_matrix
dist_sparse = csr_matrix((dist,(row,col)), shape=(num_node,num_node))

# calculate the shortest path
from scipy.sparse.csgraph import dijkstra
dist, pred = dijkstra(dist_sparse, indices = start_node, return_predecessors=True)
print(f'distance from node #{start_node:0d} to node #{end_node:0d}: {dist[end_node]:4.2f}')

# obtain the shortest path
path = []
i=end_node
if np.isinf(dist[end_node]):
    print('the path does not exist!')
else:
    while i!=start_node:
        path.append(i)
        i = pred[i]
    path.append(start_node)
    print('path=',path[::-1])

opt_path = np.asarray(path[::-1])

# plot all paths, obstacles
plt.figure(1)
#plt.triplot(xy_points[:,0], xy_points[:,1], tri.simplices)
plt.plot(x_obs_0,y_obs_0,'r',linewidth=4)
plt.plot(x_obs_1,y_obs_1,'r',linewidth=4)
plt.plot(xy_start[0],xy_start[1],'x')
plt.plot(xy_dest[0],xy_dest[1],'o')
plt.text(0,0.2,'Initial Location')
plt.text(7.5,4.2,'Destination')
plt.axis([-0.2,10,-0.2,5])
#plt.savefig('../book_v0_3/figures/python_nonconvex_obstacle_path_planning.pdf',dpi=250)

# plot the optimal path
for idx in range(0,opt_path.size-1,1):
    p_idx = opt_path[idx]
    q_idx = opt_path[idx+1]
    op_xx = [xy_points[p_idx,0], xy_points[q_idx,0]]
    op_yy = [xy_points[p_idx,1], xy_points[q_idx,1]]
    plt.plot(op_xx,op_yy,'g-.',linewidth=2)
    
#------------------------------------------------------------------------------
# resampling
#------------------------------------------------------------------------------
xy_opt_points = xy_points[opt_path,:]
dxy_opt_dist = np.sqrt(np.sum((xy_opt_points[0:-1]-xy_opt_points[1::])**2,1))
N_new_samp = 1000

xy_samp = np.empty((0,2))

for crd, dst in zip(xy_opt_points,dxy_opt_dist):
    xy_samp = np.append(xy_samp,crd + np.random.randn(N_new_samp,2)*dst,axis=0)

# found points are not inside the circular obstacle
mask_0 = ~Obs_0.contains_points(xy_samp)
xy_samp = xy_samp[mask_0,::]

mask_1 = ~Obs_1.contains_points(xy_samp)
xy_samp = xy_samp[mask_1,::]

# add start and destination node at 0 and 1st row
xy_samp = np.vstack((xy_start,xy_dest,xy_samp))

# construct graph
tri_updated = Delaunay(xy_samp)

# found triangle definition index
temp_idx=tri_updated.simplices[::,0]
temp_jdx=tri_updated.simplices[::,1]
temp_kdx=tri_updated.simplices[::,2]

# calculate the distance
dist_ij = np.sqrt(np.sum((xy_samp[temp_idx,::]-xy_samp[temp_jdx,::])**2,1))
dist_jk = np.sqrt(np.sum((xy_samp[temp_jdx,::]-xy_samp[temp_kdx,::])**2,1))
dist_ki = np.sqrt(np.sum((xy_samp[temp_kdx,::]-xy_samp[temp_idx,::])**2,1))
cut_dist = 0.1*np.sqrt(np.max(np.sum((xy_points[opt_path[0:-1],:]-xy_points[opt_path[1::],:])**2,axis=1)))

# distance thereshold for removing longer paths
cut_mask_ij = dist_ij<cut_dist
cut_mask_jk = dist_jk<cut_dist
cut_mask_ki = dist_ki<cut_dist
temp_xy_ij = np.vstack((temp_idx[cut_mask_ij],temp_jdx[cut_mask_ij]))
temp_xy_jk = np.vstack((temp_jdx[cut_mask_jk],temp_kdx[cut_mask_jk]))
temp_xy_ki = np.vstack((temp_kdx[cut_mask_ki],temp_idx[cut_mask_ki]))

# corresponding distance to the paths
dist_ij = dist_ij[cut_mask_ij]
dist_jk = dist_jk[cut_mask_jk]
dist_ki = dist_ki[cut_mask_ki]

# change format into row, column and the distance
xy_index = np.hstack((temp_xy_ij,temp_xy_jk,temp_xy_ki)).T
row_org = xy_index[::,0]
col_org = xy_index[::,1]
row = np.hstack((row_org,col_org))
col = np.hstack((col_org,row_org))
dist = np.hstack((dist_ij,dist_jk,dist_ki))
dist = np.hstack((dist,dist))
num_node = xy_samp.shape[0]

dist_sparse = csr_matrix((dist,(row,col)), shape=(num_node,num_node))
dist, pred = dijkstra(dist_sparse, indices = start_node, return_predecessors=True)

# obtain the shortest path
path = []
i=end_node
if np.isinf(dist[end_node]):
    print('the path does not exist!')
else:
    while i!=start_node:
        path.append(i)
        i = pred[i]
    path.append(start_node)
    print('path=',path[::-1])

opt_path_updated = np.asarray(path[::-1])

# plot all paths, obstacles
plt.figure(2)
plt.triplot(xy_samp[:,0], xy_samp[:,1], tri_updated.simplices)
plt.plot(x_obs_0,y_obs_0,'r',linewidth=4)
plt.plot(x_obs_1,y_obs_1,'r',linewidth=4)
plt.plot(xy_start[0],xy_start[1],'x')
plt.plot(xy_dest[0],xy_dest[1],'o')
plt.text(0,0.2,'Initial Location')
plt.text(7.5,4.2,'Destination')
plt.axis([-0.2,10,-0.2,5])

# plot the updated optimal path
plt.figure(1)
for idx in range(0,opt_path_updated.size-1,1):
    p_idx = opt_path_updated[idx]
    q_idx = opt_path_updated[idx+1]
    op_xx = [xy_samp[p_idx,0], xy_samp[q_idx,0]]
    op_yy = [xy_samp[p_idx,1], xy_samp[q_idx,1]]
    plt.plot(op_xx,op_yy,'m',linewidth=4)
    
plt.savefig('../book_v0_3/figures/python_org_resample_path_planning.pdf',dpi=250)

plt.figure(2)
for idx in range(0,opt_path_updated.size-1,1):
    p_idx = opt_path_updated[idx]
    q_idx = opt_path_updated[idx+1]
    op_xx = [xy_samp[p_idx,0], xy_samp[q_idx,0]]
    op_yy = [xy_samp[p_idx,1], xy_samp[q_idx,1]]
    plt.plot(op_xx,op_yy,'m',linewidth=4)

plt.savefig('../book_v0_3/figures/python_resample_path_planning.pdf',dpi=250)
