clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
%% Teste verschiedene Materialparameter
maxOrder = 2;
hVec = 1./(2.^(5:6)); % 32, 64
nuVec = [0.45, 0.49];%[0.05, 0.2, 0.4, 0.45, 0.49];
nComparisons = maxOrder * length(hVec) * length(nuVec);
figGrid = figure("Name","Triangulierung",'NumberTitle','off');
tiledlayout(1,nComparisons);
figSolution = figure("Name","Loesung des Elastizitaetproblems",'NumberTitle','off');
tiledlayout(2,nComparisons,TileIndexing='columnmajor');
figDefVec = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
tiledlayout(1,nComparisons);
figDefPol = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
tiledlayout(2,nComparisons,TileIndexing='columnmajor');
%TODO plotfunktionen mehr infos mitgeben

for i = 1:length(nuVec)
    nu = nuVec(i);
    for j = 1:length(hVec)
        h = hVec(j);
        for order = 1:maxOrder
            [vert,tri] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen
            [vert,tri] = extendGridLagr(vert,tri,order);
            dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
            grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
            plotGridDirichlet(grid,[],figGrid,"Triangulierung der Ordnung 1");

            % PDE
            E = 210;  % Materialparameter
            f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
            gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
            
            [U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

            deformed_area = vert; % Deformierte Liste initialisieren
            deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
            deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung
            
            % Plots
            figSolution = plotVectorfieldSolution(vert,tri,U,V,[],figSolution);
            figDefVec = plotDeformationVectors(order,vert,deformed_area,U,V,[],figDefVec);
            figDefPol = plotDeformationPolygons(tri,vert,order,deformed_area,[],figDefPol);
        end
    end
end