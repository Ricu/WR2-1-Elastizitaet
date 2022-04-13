clear; clc;
addpath('libs')
addpath('libs/distmesh')
%%
[vert,tri] = genMeshSquare(1,16);
dirichlet = (vert(:,1) == 0);
figure()
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r")
legend("Triangulierung","Dirichletrand Knoten")
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet);

E = 210; 
nu = 0.3;
f = @(vert,y) [ones(size(vert));ones(size(y))];
gD = @(x) 0*x;

[U,V] = elastSolver(grid,E,nu,f,gD);

figure()
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1")
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2")