clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
%% Mesh
h = 1/16;
[vert1,tri1] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen
[vert2,tri2] = extendGridLagr(vert1,tri1,2);
dirichlet1 = (vert1(:,1) == 0); % Dirichletrand, logischer Vektor
dirichlet2 = (vert2(:,1) == 0); % Dirichletrand, logischer Vektor
grid = struct("vert",vert1,"tri",tri1,"dirichlet",dirichlet1); % Gitter in eine Structure  bringen. 
grid2 = struct("vert",vert2,"tri",tri2,"dirichlet",dirichlet2); % Gitter in eine Structure  bringen. 

figGrid = plotGridDirichlet(grid,2,[],"Triangulierung der Ordnung 1");
figGrid = plotGridDirichlet(grid2,[],figGrid,"Triangulierung der Ordnung 2");
%% PDE aus Teil 1
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) 100*[ones(size(x));zeros(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]


order=1;    %Grad der Basisfunktionen festlegen
[U1,V1] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
order=2;    %Grad der Basisfunktionen festlegen
[U2,V2] = elastSolver(grid2,E,nu,f,gD,order); % Problem loesen


%% Loesung plotten
figSolutionPart1 = plotVectorfieldSolution(vert1,tri1,U1,V1,2);
figSolutionPart1 = plotVectorfieldSolution(vert2,tri2,U2,V2,2,figSolutionPart1);

%% Deformierte Flaeche darstellen
deformed_area1 = vert1; % Deformierte Liste initialisieren
deformed_area1(:,1) = deformed_area1(:,1) + U1; % Deformierung in x_1 Richtung
deformed_area1(:,2) = deformed_area1(:,2) + V1; % Deformierung in x_2 Richtung

figDefVecPart1 = plotDeformationVectors(1,vert1,deformed_area1,U1,V1,2);
figDefPolPart1 = plotDeformationPolygons(tri1,vert1,1,deformed_area1,2);

deformed_area2 = vert2; % Deformierte Liste initialisieren
deformed_area2(:,1) = deformed_area2(:,1) + U2; % Deformierung in x_1 Richtung
deformed_area2(:,2) = deformed_area2(:,2) + V2; % Deformierung in x_2 Richtung

figDefVecPart1 = plotDeformationVectors(2,vert2,deformed_area2,U2,V2,2,figDefVecPart1);
figDefPolPart1 = plotDeformationPolygons(tri2,vert2,2,deformed_area2,2,figDefPolPart1);
