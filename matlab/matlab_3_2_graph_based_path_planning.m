% MIT License
% 
% Copyright (c) 2022 Jongrae.K
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear;

% sparse matrix from full matrix
A_path_graph_full = [0 2 0 5 0 0;
                     2 0 1 3 4 0;
                     0 1 0 4 0 0;
                     5 3 4 0 2 4;
                     0 4 0 2 0 1;
                     0 0 0 4 1 0];
                 
                
A_path_graph_sparse = sparse(A_path_graph_full);

% sparse matrix from (row,column,values)
row_size = 6;
col_size = 6;

row = [1 1 2 2 2 3 4 4 5];
col = [2 4 3 4 5 4 5 6 6];
val = [2 5 1 3 4 4 2 4 1];

A_path_graph_sparse = sparse(row,col,val,row_size,col_size);
A_path_graph_sparse =  A_path_graph_sparse + A_path_graph_sparse';

% Shortest path planning
st_node = 1;
end_node = 6;

G_path_graph = graph(row,col,val);
[opt_path,opt_dst] = shortestpath(G_path_graph,st_node,end_node);

% Plot the result
G_graph_plot = plot(G_path_graph, ...
        'EdgeLabel',G_path_graph.Edges.Weight, ...
        'NodeFontSize',14,'EdgeFontSize',12);
highlight(G_graph_plot,opt_path,'EdgeColor','r', ...
        'LineWidth',2,'LineStyle','--');
