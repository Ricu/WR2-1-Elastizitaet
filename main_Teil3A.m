clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

%% Maximalen Grad der Basisfunktionen festlegen
maxOrder=2; 

%% Gitter erstellen (Vorarbeit)
h = 1/16; % Gitterfeinheit
[vert,tri] = genMeshSquare(1,1/h); % Knotenliste und Elementeliste erstellen

%% PDE
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
% E = 21*10^5; nu = 0.28; % Stahl
% E = 0.037; nu = 0.485; % Radiergummi

%% Figures
% Macht die Uebergabe einfacher und dient als logische Einheit
figGrid = plotGridDirichlet([],2);
figSolution = plotVectorfieldSolution([],[],[],[],2);
figDefVec = plotDeformationVectors([],[],[],[],2);
figDefPol = plotDeformationPolygons([],[],[],[],2);

for order = 1:maxOrder
    %% Restliches Gitter erstellen und und Triangulierung plotten
    [vert,tri] = extendGridLagr(vert,tri,order); % Fuer hoehere Ordnung als P1: Hinzufuegen von Knoten
    dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
    grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
    figGrid = plotGridDirichlet(grid,2,figGrid,sprintf("Triangulierung der Ordnung %i",order));

    %% Problem loesen
    [U,V] = elastSolver(grid,E,nu,f,gD,order);
    
    %% Loesung plotten
    figSolution = plotVectorfieldSolution(vert,tri,U,V,2,figSolution);

    %% Deformierte Flaeche darstellen
    % Plot der veraenderten Position der Knoten mit Verschiebungsvektoren
    [figDefVec,~,defVert] = plotDeformationVectors(vert,U,V,order,2,figDefVec);
    % Plot des defomierten Polygons
    figDefPol = plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol,sprintf("Ordnung %g:",order));
end
