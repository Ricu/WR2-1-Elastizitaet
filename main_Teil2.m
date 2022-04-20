clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Polygon 1
pv = [-0.4,-0.5;...
       0.4,-0.2;...
       0.4,-0.7;...
       1.5,-0.4;...
       0.9, 0.1;...
       1.6, 0.8;...
       0.5, 0.5;...
       0.2, 1  ;...
       0.1, 0.4;...
      -0.4, 0.7;... % Anders als im beispiel
      -0.4,-0.5]; % Eckpunkte definieren
[vert,tri] = distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],true,pv,pv); % Gitter erzeugen

%% Polygon 2
pv = [ 0,  0;...
       0, 22;...
      24, 30;...
      24, 22;...
       0,  0]; % Eckpunkte definieren
distmesh2d(@dpoly,@huniform,6,[0,0; 24,30],true,pv,pv); % Gitter erzeugen %Verusacht endlosschleife

%% Dirichletknoten hinzufuegen und plotten
dirichlet = (vert(:,1) == -0.4); % Dirichletrand, logischer Vektor
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit
plotGridDirichlet(grid,1,[],"Triangulierung der Ordnung 1");

%% Testen
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion
order = 1;    %Grad der Basisfunktionen festlegen

[U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

figSolution = plotVectorfieldSolution(vert,tri,U,V,1);

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung

plotDeformationVectors(order,vert,deformed_area,U,V,1);
plotDeformationPolygons(tri,vert,order,deformed_area,1);