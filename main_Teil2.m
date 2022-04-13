clear; clc;
addpath('libs')
addpath('libs/distmesh')

%% Polygon 1
pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
       1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
[vert,tri] = distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);

%% Polygon 2
% pv = [0, 0; 0, 22; 24, 30; 24, 22];
% distmesh2d(@dpoly,@huniform,6,[0,0; 24,30],pv,pv); %Verusacht endlosschleife

%% Dirichletknoten hinzufügen und plotten
dirichlet = (vert(:,1) == -0.4);
figure()
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r")
legend("Triangulierung","Dirichletrand Knoten")
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet);

%% Testen
E = 210; 
nu = 0.3;
f = @(vert,y) [ones(size(vert));ones(size(y))];
gD = @(x) 0*x;

[U,V] = elastSolver(grid,E,nu,f,gD);

figure()
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1")
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2")