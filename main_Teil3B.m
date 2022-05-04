clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

%% Verschiedene Werte der Querkontraktionszahl testen (nu->1/2)
nuVec = [0.45, ...
    0.49999999999999, ...
    0.499999999999999, ...
    0.4999999999999999];
% Ergebnisse:
% Komische Triangulierung in P2 für 0.49999999999999
% Komische Triangulierung in P1 für 0.499999999999999
% Chaos fuer 0.4999999999999999
% Sehr kaputte Triangulierung für 0.49999999999999999

%% Maximalen Grad der Basisfunktionen festlegen
maxOrder=2;

%% PDE
E = 40;  % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]

%% Figures
figDefPol = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
tiledlayout(length(nuVec),2,'TileSpacing','Compact','Padding','Compact',TileIndexing='rowmajor');
axes = cell(length(nuVec),2);

%% Loop
for i = 1:length(nuVec)
    nu = nuVec(i);
    for order = 1:maxOrder
        %% Gitter erstellen und und Triangulierung plotten
        h = 1/16; % Gitterfeinheit
        [vert,tri] = genMeshSquare(1,1/h); % Knotenliste und Elementeliste erstellen
        [vert,tri] = extendGridLagr(vert,tri,order); % Fuer hoehere Ordnung als P1: Hinzufuegen von Knoten
        dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen
        
        %% Problem loesen
        [U,V] = elastSolver(grid,E,nu,f,gD,order);
        
        %% Deformierte Flaeche darstellen
        defVert = vert; % Deformierte Liste initialisieren
        defVert(:,1) = defVert(:,1) + U; % Deformierung in x_1 Richtung
        defVert(:,2) = defVert(:,2) + V; % Deformierung in x_2 Richtung
        % Plot des defomierten Polygons
        customTitle = sprintf("nu = %.17f, order = %i",nu,order);
        plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol,customTitle);
        axes{i,order} = gca;
    end
end