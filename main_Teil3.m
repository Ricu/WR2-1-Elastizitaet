clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
%% Mesh
h = 1/16;
[vert1,tri1] = genMeshSquare(1,1/h); % TODO val Punkte und Dreiecke erstellen
[vert2,tri2] = extendGridLagr(vert1,tri1,2);
dirichlet1 = (vert1(:,1) == 0); % Dirichletrand, logischer Vektor
dirichlet2 = (vert2(:,1) == 0);
