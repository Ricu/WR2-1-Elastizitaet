clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Mesh
h = 1/16;
[vert,tri] = genMeshSquare(1,1/h); % Knoten und Elemente erstellen
dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit
plotGridDirichlet(grid,1,[],"Triangulierung der Ordnung 1");

%% PDE 
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
order = 1;    %Grad der Basisfunktionen festlegen

[U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

%% Loesung plotten
figSolution = plotVectorfieldSolution(vert,tri,U,V,1);

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung

[~,~,defVert] = plotDeformationVectors(vert,U,V,order,1);
plotDeformationPolygons(vert,tri,deformed_area,order,1);