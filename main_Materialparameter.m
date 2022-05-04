clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Teste verschiedene Materialparameter
% Dieses Skript dient dazu, die durchschnittliche Aenderung der
% Deformation in den Knoten in Abhaengigkeit von der Schrittweite bzw dem
% Materialparameter nu zu analysieren
maxOrder = 1;
hVec = 1./(2.^(2:6));  % 4, 8, 16, 32
nuVec = [0.1, 0.3, 0.45, 0.49];

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
            verts{i} = vert;

            % PDE
            E = 210;  % Materialparameter
            f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
            gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
            
            [U{i,j},V{i,j}] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen
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