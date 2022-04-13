clear; clc; % Konsolen Output und Variablen löschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Polygon 1
pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
       1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]; % Eckpunkte definieren
[vert,tri] = distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv); % Gitter erzeugen

%% Polygon 2
% pv = [0, 0; 0, 22; 24, 30; 24, 22]; % Eckpunkte definieren
% distmesh2d(@dpoly,@huniform,6,[0,0; 24,30],pv,pv); % Gitter erzeugen %Verusacht endlosschleife

%% Dirichletknoten hinzufügen und plotten
dirichlet = (vert(:,1) == -0.4); % Dirichletrand, logischer Vektor
figure() % Neues Fenster erzeugen
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufügen
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die Übergabe einfacher und dient als logische Einheit

%% Testen
E = 210; nu = 0.3; % Materialparameter
f = @(vert,y) [ones(size(vert));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion

[U,V] = elastSolver(grid,E,nu,f,gD); % Problem lösen

figure() % Neues Fenster erzeugen
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1") % Lösung in x_1 Richtung plotten
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2") % Lösung in x_2 Richtung plotten