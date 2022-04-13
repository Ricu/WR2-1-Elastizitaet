clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%%
[vert,tri] = genMeshSquare(1,16); % Punkte und Dreiecke erstellen
dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
figure() % Neues Fenster erzeugen
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufuegen
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit

E = 210; nu = 0.3; % Materialparameter
f = @(vert,y) [ones(size(vert));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion

% func = @() elastSolver(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Optimiert",h, timeit(func))
% 
% func = @() elastSolver2(grid,E,nu,f,gD);
% fprintf("%15s: Benoetigte Zeit fuer 1/h = %i: %fs\n", "Nicht Optimiert",h, timeit(func))

[U,V] = elastSolver(grid,E,nu,f,gD); % Problem loesen

%% Loesung plotten
figure() % Neues Fenster erzeugen
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1") % Loesung in x_1 Richtung plotten
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2") % Loesung in x_2 Richtung plotten

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung

figure() % Neues Fenster erzeugen
scatter(vert(:,1),vert(:,2),'filled'); hold on; % Urspruengliche Flaeche plotten
scatter(deformed_area(:,1),deformed_area(:,2)); % Deformierte Flaeche plotten
legend("Old","New") % Legende hinzufügen
quiver(vert(:,1),vert(:,2),U,V,0) % Berechnetes Vektorfeld (Verschiebung) plotten
