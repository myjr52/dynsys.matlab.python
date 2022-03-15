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

# recycle delaunay code by chaining to voronoi
# construct graph using voronoi
from scipy.spatial import Voronoi
vor = Voronoi(xy_points)
xy_points = vor.vertices
num_sample = xy_points.shape[0]
xy_points_org_4_draw = xy_points

#-----------------------------------------------------------------------------
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
mask_1 = ~Obs_1.contains_points(xy_points)
mask_all = np.logical_and(mask_0,mask_1)
#-----------------------------------------------------------------------------

# edges
temp_index = np.array(vor.ridge_vertices)

# remove edges inside the boundary: this part needs to be updated for
# faster execusion
temp_index_update = []
set_index = np.arange(num_sample)
set_index = set_index[~mask_all]
for add_idx_p1_p2 in enumerate(temp_index):
    add_idx = add_idx_p1_p2[0]
    p1_idx = add_idx_p1_p2[1][0]
    p2_idx = add_idx_p1_p2[1][1]
    if not(np.isin(p1_idx,set_index) or np.isin(p2_idx,set_index)):
        temp_index_update.append(temp_index[add_idx])
#--------------------------------------------------------------
        
temp_index_update = np.array(temp_index_update)
temp_idx = temp_index_update[:,0]
temp_jdx = temp_index_update[:,1]

# remove edges outside of the boundary
cut_mask = temp_idx <= num_sample-1
temp_idx = temp_idx[cut_mask]
temp_jdx = temp_jdx[cut_mask]
cut_mask = temp_jdx <= num_sample-1
temp_idx = temp_idx[cut_mask]
temp_jdx = temp_jdx[cut_mask]

cut_mask = temp_idx >= 0
temp_idx = temp_idx[cut_mask]
temp_jdx = temp_jdx[cut_mask]
cut_mask = temp_jdx >= 0
temp_idx = temp_idx[cut_mask]
temp_jdx = temp_jdx[cut_mask]    

dd_all = np.sqrt(np.sum((xy_points[temp_idx]-xy_points[temp_jdx])**2,1))
cut_dist = np.mean(dd_all)+1*np.std(dd_all)

# distance thereshold for removing longer paths
cut_mask_ij = dd_all<cut_dist
temp_xy_ij = np.vstack((temp_idx[cut_mask_ij],temp_jdx[cut_mask_ij]))

# corresponding distance to the paths
dist_ij = dd_all[cut_mask_ij]

# change format into row, column and the distance
xy_index = temp_xy_ij.T
row_org = xy_index[:,0]
col_org = xy_index[:,1]
row = np.hstack((row_org,col_org))
col = np.hstack((col_org,row_org))
dist = dist_ij
dist = np.hstack((dist,dist))
num_node = xy_points.shape[0]

# construct the distance matrix
from scipy.sparse import csr_matrix
dist_sparse = csr_matrix((dist,(row,col)), shape=(num_node,num_node))

# update start_node and end_node to the closest points in the list
start_node = np.argmin(np.sum((xy_points-xy_start)**2,1))
end_node = np.argmin(np.sum((xy_points-xy_dest)**2,1))


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
plt.figure()
x_p1_p2=np.array([xy_points_org_4_draw[temp_xy_ij[0]][:,0], xy_points_org_4_draw[temp_xy_ij[1]][:,0]])
y_p1_p2=np.array([xy_points_org_4_draw[temp_xy_ij[0]][:,1], xy_points_org_4_draw[temp_xy_ij[1]][:,1]])
plt.plot(x_p1_p2,y_p1_p2,'b-')    
plt.plot(x_obs_0,y_obs_0,'r',linewidth=4)
plt.plot(x_obs_1,y_obs_1,'r',linewidth=4)
plt.plot(xy_start[0],xy_start[1],'x')
plt.plot(xy_dest[0],xy_dest[1],'o')
plt.text(0,0.2,'Initial Location')
plt.text(7.5,4.2,'Destination')
plt.axis([-0.2,10,-0.2,5])

# plot the optimal path
for idx in range(0,opt_path.size-1,1):
    p_idx = opt_path[idx]
    q_idx = opt_path[idx+1]
    op_xx = [xy_points[p_idx,0], xy_points[q_idx,0]]
    op_yy = [xy_points[p_idx,1], xy_points[q_idx,1]]
    plt.plot(op_xx,op_yy,'g',linewidth=4)
    
# connect the starting point and the end point to the path
op_xx = [xy_points[start_node,0], xy_start[0]]
op_yy = [xy_points[start_node,1], xy_start[1]]
plt.plot(op_xx,op_yy,'g',linewidth=4)
op_xx = [xy_points[end_node,0], xy_dest[0]]
op_yy = [xy_points[end_node,1], xy_dest[1]]
plt.plot(op_xx,op_yy,'g',linewidth=4)
