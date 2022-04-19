clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
%% Mesh
h = 1/16;
[vert1,tri1] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen
[vert2,tri2] = extendGridLagr(vert1,tri1,2);
dirichlet1 = (vert1(:,1) == 0); % Dirichletrand, logischer Vektor
dirichlet2 = (vert2(:,1) == 0);

figGrid = plotGridDirichlet(vert1,tri1,dirichlet1,2);
figGrid = plotGridDirichlet(vert2,tri2,dirichlet2,2,figGrid);

grid = struct("vert",vert1,"tri",tri1,"dirichlet",dirichlet1); % Gitter in eine Structure  bringen. 
grid2 = struct("vert",vert2,"tri",tri2,"dirichlet",dirichlet2); % Gitter in eine Structure  bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit

%% PDE aus Teil 1
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]


% func = @() elastSolver(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Optimiert",h, timeit(func))
% 
% func = @() elastSolver2(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Nicht Optimiert",h, timeit(func))

order=1;    %Grad der Basisfunktionen festlegen
[U1,V1] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
order=2;    %Grad der Basisfunktionen festlegen
% [U2,V2] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
U = U1;
V = V1;

%% Loesung plotten
figSolution = plotVectorfieldSolution(vert1,tri1,U1,V1,nComparisons);
figSolution = plotVectorfieldSolution(vert2,tri2,U2,V2,nComparisons,figSolution);

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U1; % Deformierung in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V1; % Deformierung in x_2 Richtung

plotDeformationVectors(vert1,deformed_area,U1,V1,1)

plotDeformationPolygons(vert1,deformed_area,1)