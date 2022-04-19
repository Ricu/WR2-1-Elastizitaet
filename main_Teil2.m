clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Polygon 1
pv = [-0.4,-0.5;...
       0.4,-0.2;...
       0.4,-0.7;...
       1.5,-0.4;...
       0.9, 0.1;...
       1.6, 0.8;...
       0.5, 0.5;...
       0.2, 1  ;...
       0.1, 0.4;...
      -0.4, 0.7;... % Anders als im beispiel
      -0.4,-0.5]; % Eckpunkte definieren
[vert,tri] = distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],true,pv,pv); % Gitter erzeugen

%% Polygon 2
pv = [ 0,  0;...
       0, 22;...
      24, 30;...
      24, 22;...
       0,  0]; % Eckpunkte definieren
distmesh2d(@dpoly,@huniform,6,[0,0; 24,30],true,pv,pv); % Gitter erzeugen %Verusacht endlosschleife

%% Dirichletknoten hinzufuegen und plotten
dirichlet = (vert(:,1) == -0.4); % Dirichletrand, logischer Vektor
figure('Name','Triangulierung') % Neues Fenster erzeugen
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufuegen
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit

%% Testen
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion
order=1;    %Grad der Basisfunktionen festlegen

[U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

figure('Name','Deformierung in x_1 und x_2 Richtung') % Neues Fenster erzeugen
subplot(1,2,1), trisurf(tri,vert(:,1),vert(:,2),U), title("(u_h)_1: x_1 Richtung") % Loesung in x_1 Richtung plotten
subplot(1,2,2), trisurf(tri,vert(:,1),vert(:,2),V), title("(u_h)_2: x_2 Richtung") % Loesung in x_2 Richtung plotten

%% Deformierte Flaeche darstellen
deformed_area = vert; % Deformierte Liste initialisieren
deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung

figure('Name','Gebietsvergleich: vor und nach Deformierung') % Neues Fenster erzeugen
scatter(vert(:,1),vert(:,2),'k','filled'); hold on; % Urspruengliche Flaeche plotten
scatter(deformed_area(:,1),deformed_area(:,2),46); % Deformierte Flaeche plotten
quiver(vert(:,1),vert(:,2),U,V,0) % Berechnetes Vektorfeld (Verschiebung) plotten
legend("Original","Deformiert","Verschiebung") % Legende hinzufuegen

figure('Name','Gebietsvergleich: vor und nach Deformierung')
subplot(1,2,1); patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
title('Urspruengliches Gebiet')
axis equal tight; title("Urspruengliches Gebiet")
subplot(1,2,2); patch('vertices',deformed_area,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
title('Deformiertes Gebiet')
axis equal tight; title("Deformiertes Gebiet")