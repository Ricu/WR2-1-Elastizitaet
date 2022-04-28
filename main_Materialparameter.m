clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Teste verschiedene Materialparameter
maxOrder = 1;
hVec = 1./(2.^(2:6)); % 32, 64
nuVec = [0.1, 0.3, 0.45, 0.49]; % [0.05, 0.2, 0.4, 0.45, 0.49];
% nComparisons = maxOrder * length(hVec) * length(nuVec);
% xN = length(nuVec);
% yN = length(hVec)*maxOrder;
% % figGrid = figure("Name","Triangulierung",'NumberTitle','off','PaperType','A4');
% % tiledlayout(xN,yN);
% figSolution = figure("Name","Loesung des Elastizitaetproblems",'NumberTitle','off');
% tiledlayout(xN,2*yN,TileIndexing='columnmajor');
% % figDefVec = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
% % tiledlayout(xN,yN);
% figDefPol = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
% tiledlayout(xN,2*yN,TileIndexing='columnmajor');
%TODO plotfunktionen mehr infos mitgeben


verts = cell(length(hVec),1);
U = cell(length(hVec),length(nuVec));
V = cell(length(hVec),length(nuVec));
diffH = cell(length(hVec),length(hVec));
diffNu = cell(length(nuVec),length(nuVec));

for i = 1:length(hVec)
    h = hVec(i);
    for j = 1:length(nuVec)
        nu = nuVec(j);
        for order = 1:maxOrder
            [vert,tri] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen
            [vert,tri] = extendGridLagr(vert,tri,order);
            dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
            grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen. 
%             plotGridDirichlet(grid,[],figGrid,"Triangulierung der Ordnung 1");
            verts{i} = vert;

            % PDE
            E = 210;  % Materialparameter
            f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
            gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
            
            [U{i,j},V{i,j}] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

%             deformed_area = vert; % Deformierte Liste initialisieren
%             deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
%             deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung
            
%             % Plots
%             figSolution = plotVectorfieldSolution(vert,tri,U,V,[],figSolution);
% %             figDefVec = plotDeformationVectors(order,vert,deformed_area,U,V,[],figDefVec);
%             figDefPol = plotDeformationPolygons(tri,vert,order,deformed_area,[],figDefPol);

%                         figSolution = plotVectorfieldSolution(vert,tri,U,V,1);
%             figDefVec = plotDeformationVectors(order,vert,deformed_area,U,V,1);
%             figDefPol = plotDeformationPolygons(tri,vert,order,deformed_area,1);

            
        end
    end
end

for i = 1: length(hVec)
    for j = 1:i-1
        [ind,loc] = ismember(verts{i},verts{j},'rows');
        % index j: weniger feines gitter
        diffH{i,j} = [sum(U{i,3}(ind) - U{j,3})/numel(U{j,3}),sum(V{i,3}(ind) - V{j,3})/numel(V{j,3})];
    end
end

for i = 1: length(nuVec)
    for j = 1:i-1
        [ind,loc] = ismember(verts{3},verts{3},'rows');
        % index j: weniger feines gitter
        diffNu{i,j} = [sum(U{3,i}(ind) - U{3,j})/numel(U{3,j}),sum(V{3,i}(ind) - V{3,j})/numel(V{3,j})];
    end
end