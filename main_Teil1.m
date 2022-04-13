clear; clc; % Konsolen Output und Variablen löschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%%
[vert,tri] = genMeshSquare(1,16); % Punkte und Dreiecke erstellen
dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
figure() % Neues Fenster erzeugen
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufügen
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die Übergabe einfacher und dient als logische Einheit

E = 210; nu = 0.3; % Materialparameter
f = @(vert,y) [ones(size(vert));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion

% func = @() elastSolver(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Optimiert",h, timeit(func))
% 
% func = @() elastSolver2(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Nicht Optimiert",h, timeit(func))

[U,V] = elastSolver(grid,E,nu,f,gD); % Problem lösen

figure() % Neues Fenster erzeugen
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1") % Lösung in x_1 Richtung plotten
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2") % Lösung in x_2 Richtung plotten