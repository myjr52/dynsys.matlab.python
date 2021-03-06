#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Jongrae.K

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
from scipy import sparse

# sparse matrix from full matrix
A_path_graph_full = np.array([[0, 2, 0, 5, 0, 0],
                              [2, 0, 1, 3, 4, 0],
                              [0, 1, 0, 4, 0, 0],
                              [5, 3, 4, 0, 2, 4],
                              [0, 4, 0, 2, 0, 1],
                              [0, 0, 0, 4, 1, 0]])
                 
A_path_graph_sparse = sparse.csc_matrix(A_path_graph_full)


# sparse matrix from (row,column,values)
row_size = 6
col_size = 6

row = np.array([0, 0, 1, 1, 1, 2, 3, 3, 4])
col = np.array([1, 3, 2, 3, 4, 3, 4, 5, 5])
val = np.array([2, 5, 1, 3, 4, 4, 2, 4, 1])

A_path_graph_sparse = sparse.csc_matrix((val,(row,col)),shape=(row_size,col_size))
A_path_graph_sparse = A_path_graph_sparse+A_path_graph_sparse.transpose()

# shortest path from node 1
start_node = 0
end_node = 5

from scipy.sparse.csgraph import dijkstra
dist, pred = dijkstra(A_path_graph_sparse, indices = start_node, return_predecessors=True)

# print out the distance from start_node to end_node
print(f"distance from {start_node} to {end_node}: {dist[end_node]}.")

# construct the path
path = []
idx=end_node
while idx!=start_node:
    path.append(idx)
    idx = pred[idx]
    
path.append(start_node)
path.reverse()
print('path=',path)
opt_path_edge = [(i,j) for i,j in zip(path[0:-1],path[1::])]

# draw graph
import matplotlib.pyplot as plt
fig, ax = plt.subplots(nrows=1,ncols=1)

import networkx as nx
A_graph_nx = nx.from_scipy_sparse_matrix(A_path_graph_sparse)
pos=nx.spring_layout(A_graph_nx)
pos=nx.planar_layout(A_graph_nx)
edge_labels=nx.get_edge_attributes(A_graph_nx,'weight')
edge_color = ['red' if key in opt_path_edge else 'green' for key in edge_labels.keys()]
nx.draw(A_graph_nx, pos, node_size=500, node_color='yellow', edge_color=edge_color, labels={node:node for node in A_graph_nx.nodes()})
nx.draw_networkx_edge_labels(A_graph_nx,pos=pos,edge_labels=edge_labels)

#fig.set_size_inches(9,6)    
#fig.savefig('../book_v0_3/figures/python_graph_shortest_path_planning.pdf',dpi=250)
