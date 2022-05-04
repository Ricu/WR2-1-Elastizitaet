clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

%% Gitter erstellen und und Triangulierung plotten
h = 1/16; % Gitterfeinheit
[vert,tri] = genMeshSquare(1,1/h); % Knotenliste und Elementeliste erstellen
dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
% Gitter in eine Structure  bringen
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet);  
% Macht die Uebergabe einfacher und dient als logische Einheit
figDir = plotGridDirichlet(grid,1,[],"Triangulierung der Ordnung 1");

%% PDE 
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
order = 1;    % Grad der Basisfunktionen festlegen

%% Problem loesen
[U,V] = elastSolver(grid,E,nu,f,gD,order);

%% Loesung plotten
figSolution = plotVectorfieldSolution(vert,tri,U,V,1);

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U; % Deformation in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V; % Deformation in x_2 Richtung

% Plot der veraenderten Position der Knoten mit Verschiebungsvektoren
[figDefVec,~,defVert] = plotDeformationVectors(vert,U,V,order,1);
% Plot des defomierten Polygons
figDefPol = plotDeformationPolygons(vert,tri,deformed_area,order,1); 