clear;

%%
% number of samples
num_sample = 1000;

% map size
x_min = 0; x_max = 10;
y_min = 0; y_max = 5;

% starting point
xy_start = [0 0];
xy_dest  = [9 4];

% spread num_sample random points over the map area
xn=rand(1,num_sample)*(x_max-x_min) + x_min;
yn=rand(1,num_sample)*(y_max-y_min) + y_min;

%% recycle voronoi graph by change voronoi
% divide region using voronoi
% [vx,vy] = voronoi(xn,yn); 

%% to delaunay
% divide region using delaunday
triangle_v_idx = delaunay(xn,yn);
v_idx = triangle_v_idx(:,1);
v_jdx = triangle_v_idx(:,2);
v_kdx = triangle_v_idx(:,3);

vx = [ xn(v_idx) xn(v_jdx) xn(v_kdx); 
       xn(v_jdx) xn(v_kdx) xn(v_idx)];
vy = [ yn(v_idx) yn(v_jdx) yn(v_kdx); 
       yn(v_jdx) yn(v_kdx) yn(v_idx)];

% remove long edges
dd_all = sqrt((vx(1,:)-vx(2,:)).^2 + (vy(1,:)-vy(2,:)).^2);
cut_dist = mean(dd_all) + 2*std(dd_all);
cut_idx = dd_all < cut_dist;
vx = vx(:,cut_idx);
vy = vy(:,cut_idx);
 
%%
% reject points outside the map region
idx = (vx(1,:) < x_min) | (vx(2,:) < x_min);
vx(:,idx) = [];
vy(:,idx) = [];
idx = (vx(1,:) > x_max) | (vx(2,:) > x_max);
vx(:,idx) = [];
vy(:,idx) = [];
idx = (vy(1,:) < y_min) | (vy(2,:) < y_min);
vx(:,idx) = [];
vy(:,idx) = [];
idx = (vy(1,:) > y_max) | (vy(2,:) > y_max);
vx(:,idx) = [];
vy(:,idx) = [];

% circular obstacle
th=0:0.01:2*pi;
c_cx = 3; c_cy = 3; c_r = 1.5;
xc=c_r*cos(th)+c_cx;
yc=c_r*sin(th)+c_cy;

% polygon obstacle
xv = [6; 8; 8; 5; 5; 7; 7; 6; 6]; 
yv = [1; 1; 4; 4; 3; 3; 2; 2; 1];

% draw sampling points & obstacles
figure(1); clf; 
plot(xn,yn,'k.'); 
hold on;
plot(xc,yc,'r-','LineWidth',2); 
plot(xv,yv,'r-','LineWidth',2);
axis equal;
axis([x_min-0.5 x_max y_min-0.5 y_max]);
plot(vx,vy,'b-');

%%
% remove vertices inside the circular object
vx1=vx(1,:);
vy1=vy(1,:);
r_sq = (vx1-c_cx).^2+(vy1-c_cy).^2;
idx1=(r_sq < c_r^2);

vx2=vx(2,:);
vy2=vy(2,:);
r_sq = (vx2-c_cx).^2+(vy2-c_cy).^2;
idx2=(r_sq < c_r^2);

idx= or(idx1,idx2);
vx(:,idx)=[]; 
vy(:,idx)=[];

% remove vertices inside the polygon
vx1=vx(1,:);
vy1=vy(1,:);
in = inpolygon(vx1,vy1,xv,yv);
vx(:,in)=[];
vy(:,in)=[];

vx2=vx(2,:);
vy2=vy(2,:);
in = inpolygon(vx2,vy2,xv,yv);
vx(:,in)=[];
vy(:,in)=[];

% add start & destination points to the graph
vx1d = vx(:);
vy1d = vy(:);
dr = kron(ones(length(vx1d),1),xy_start) - [vx1d vy1d];
[~,min_id] = min(sum(dr.^2,2));
vx = [vx [xy_start(1); vx1d(min_id)]];
vy = [vy [xy_start(2); vy1d(min_id)]];

dr = kron(ones(length(vx1d),1),xy_dest) - [vx1d vy1d];
[~,min_id] = min(sum(dr.^2,2));
vx = [vx [xy_dest(1); vx1d(min_id)]];
vy = [vy [xy_dest(2); vy1d(min_id)]];


% draw voronoi after removing points in the obstacles
figure(1); clf; 
plot(xn,yn,'k.'); 
hold on;
plot(xc,yc,'r-','LineWidth',2); 
plot(xv,yv,'r-','LineWidth',2);
axis equal;
plot(vx,vy,'b.-'); 
axis([x_min-0.5 x_max y_min-0.5 y_max]);
plot(xy_start(1),xy_start(2),'bx','MarkerSize',5,'LineWidth',5);
plot(xy_dest(1),xy_dest(2),'ro','MarkerSize',5,'MarkerFacecolor','red');

%% construct graph
xy_1 = [vx(1,:); vy(1,:)];
xy_2 = [vx(2,:); vy(2,:)];
xy_12 = [xy_1 xy_2]';
[node_coord,~,node_index]=unique(xy_12,'rows');
st_node_index = node_index(1:length(vx));
ed_node_index = node_index(length(vx)+1:end);
dst_edges = sqrt(sum((xy_1-xy_2).^2));
st_node = node_index(length(vx)-1);
ed_node = node_index(length(vx));

row = [st_node_index(:); ed_node_index(:)];
col = [ed_node_index(:); st_node_index(:)];
val = kron([1;1],dst_edges(:));
G_path_graph = graph(row,col,val);

%% calculate optimal path and plot the path
[opt_path_idx,opt_dst] = shortestpath(G_path_graph,st_node,ed_node);
opt_path = node_coord(opt_path_idx,:);

% draw voronoi after removing points in the obstacles
figure(1); clf; 
plot(xn,yn,'k.'); 
hold on;
plot(xc,yc,'r-','LineWidth',2); 
plot(xv,yv,'r-','LineWidth',2);
axis equal;
plot(vx,vy,'b.-'); 
axis([x_min-0.5 x_max y_min-0.5 y_max]);
plot(xy_start(1),xy_start(2),'bx','MarkerSize',5,'LineWidth',5);
plot(xy_dest(1),xy_dest(2),'ro','MarkerSize',5,'MarkerFacecolor','red');
plot(opt_path(:,1),opt_path(:,2),'g-','LineWidth',2);