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
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]


order=1;    %Grad der Basisfunktionen festlegen
[U1,V1] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
order=2;    %Grad der Basisfunktionen festlegen
[U2,V2] = elastSolver(grid2,E,nu,f,gD,order); % Problem loesen


%% Loesung plotten
figSolution = plotVectorfieldSolution(vert1,tri1,U1,V1,2);
figSolution = plotVectorfieldSolution(vert2,tri2,U2,V2,2,figSolution);

%% Deformierte Flaeche darstellen
deformed_area1 = vert1; % Deformierte Liste initialisieren
deformed_area1(:,1) = deformed_area1(:,1) + U1; % Deformierung in x_1 Richtung
deformed_area1(:,2) = deformed_area1(:,2) + V1; % Deformierung in x_2 Richtung

figDefVec = plotDeformationVectors(1,vert1,deformed_area1,U1,V1,2);
figDefPol = plotDeformationPolygons(tri1,vert1,1,deformed_area1,2);

deformed_area2 = vert2; % Deformierte Liste initialisieren
deformed_area2(:,1) = deformed_area2(:,1) + U2; % Deformierung in x_1 Richtung
deformed_area2(:,2) = deformed_area2(:,2) + V2; % Deformierung in x_2 Richtung

figDefVec = plotDeformationVectors(2,vert2,deformed_area2,U2,V2,2,figDefVec);
figDefPol = plotDeformationPolygons(tri2,vert2,2,deformed_area2,2,figDefPol);


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
            [vert,tri] = extendGridLagr(vert,tri,2);
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




