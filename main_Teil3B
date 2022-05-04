clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

h = 1/16;
maxOrder = 2;
% Komisch in P2 für 0.49999999999999
% Komisch in P1 für 0.499999999999999
% Chaos absolut für 0.4999999999999999
% Richtig kaput für 0.49999999999999999
nuVec = [0.45, ...
         0.49999999999999, ...
         0.499999999999999, ...
         0.4999999999999999];
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
        [vert,tri] = genMeshSquare(1,1/h);
        [vert,tri] = extendGridLagr(vert,tri,order);
        dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
        [U,V] = elastSolver(grid,E,nu,f,gD,order);
        
        defVert = vert; % Deformierte Liste initialisieren
        defVert(:,1) = defVert(:,1) + U; % Deformierung in x_1 Richtung
        defVert(:,2) = defVert(:,2) + V; % Deformierung in x_2 Richtung
        customTitle = sprintf("nu = %.17f, order = %i",nu,order);
        plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol,customTitle);
        axes{i,order} = gca;
    end
end