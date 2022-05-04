clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

% Dieses Skript dient dazu, die durchschnittliche Aenderung der
% Deformation in den Knoten in Abhaengigkeit von der Schrittweite bzw dem
% Materialparameter nu zu analysieren
%% Teste verschiedene Gitterfeinheiten und Werte fuer die Querkontraktionszahl
hVec = 1./(2.^(2:4));  % 1/4, 1/8, 1/16, 1/32 6
nuVec = [0.1, 0.3, 0.45, 0.49];

%% PDE
E = 210;  % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
maxOrder = 1; % Grad der Basisfunktionen festlegen

%% Loop
verts = cell(length(hVec),1);
U = cell(length(hVec),length(nuVec));
V = cell(length(hVec),length(nuVec));
for i = 1:length(hVec)
    h = hVec(i);
    for j = 1:length(nuVec)
        nu = nuVec(j);
        for order = 1:maxOrder
             %% Gitter erstelle
            [vert,tri] = genMeshSquare(1,1/h); % Knotenliste und Elementeliste erstellen
            [vert,tri] = extendGridLagr(vert,tri,order); % Fuer hoehere Ordnung als P1: Hinzufuegen von Knoten
            dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
            grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen
            verts{i} = vert;
            
             %% Problem loesen
            [U{i,j},V{i,j}] = elastSolver(grid,E,nu,f,gD,order);
        end
    end
end

%% Durchschnittliche Aenderung der Deformation in den Knoten in Abhaengigkeit von der Schrittweite berechnen
diffH = cell(length(hVec),length(hVec));
columnNames=cell(length(hVec),1);
for i = 1: length(hVec)
    for j = 1:i-1 % Index j: weniger feines Gitter
        % Pruefe, ob Knoten in  beiden Knotenlisten enthalten sind
        % ind: logischer Vektor; loc: gibt Position des Knotens in der anderen Liste an
        [ind,loc] = ismember(verts{i},verts{j},'rows');
        diffH{i,j} = [sum(U{i,3}(ind) - U{j,3})/numel(U{j,3}),sum(V{i,3}(ind) - V{j,3})/numel(V{j,3})];
    end
    columnNames{i}=sprintf('h=%g',hVec(i));
end
% Ausgabe der Ergebnisse in Tabelle 1
T_h = cell2table(diffH,"RowNames",columnNames,"VariableNames",columnNames)

%% Durchschnittliche Aenderung der Deformation in den Knoten in Abhaengigkeit des Materialparameters nu berechnen
diffNu = cell(length(nuVec),length(nuVec));
columnNames=cell(length(hVec),1);
for i = 1: length(nuVec)
    for j = 1:i-1 % Index j: kleineres nu
        % Pruefe, ob Knoten in  beiden Knotenlisten enthalten sind
        % ind: logischer Vektor; loc: gibt Position des Knotens in der anderen Liste an
        [ind,loc] = ismember(verts{3},verts{3},'rows');
        diffNu{i,j} = [sum(U{3,i}(ind) - U{3,j})/numel(U{3,j}),sum(V{3,i}(ind) - V{3,j})/numel(V{3,j})];
    end
    columnNames{i}=sprintf('nu=%g',nuVec(i));
end
% Ausgabe der Ergebnisse in Tabelle 2
T_nu = cell2table(diffNu,"RowNames",columnNames,"VariableNames",columnNames)
