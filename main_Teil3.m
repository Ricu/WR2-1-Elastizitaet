clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Mesh
h = 1/16;
[vert,tri] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen

%% PDE aus Teil 1
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
% E = 21*10^5; nu = 0.28; % Stahl
% E = 0.037; nu = 0.485; % Radiergummi

%% Figure
figGrid = plotGridDirichlet([],2);
figSolution = plotVectorfieldSolution([],[],[],[],2);
figDefVec = plotDeformationVectors([],[],[],[],2);
figDefPol = plotDeformationPolygons([],[],[],[],2);

for order = 1:2
    [vert,tri] = extendGridLagr(vert,tri,order);
    dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
    grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
    figGrid = plotGridDirichlet(grid,2,figGrid,sprintf("Triangulierung der Ordnung %i",order));

    %% Loesung plotten
    [U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
    figSolution = plotVectorfieldSolution(vert,tri,U,V,2,figSolution);

    %% Deformierte Flaeche darstellen
    [figDefVec,~,defVert] = plotDeformationVectors(vert,U,V,order,2,figDefVec);
    figDefPol = plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol);
end






