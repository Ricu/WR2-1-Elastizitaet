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
      -0.4, 0.7;...
      -0.4,-0.5]; % Eckpunkte definieren
bbox = [min(pv(:,1)), min(pv(:,2)); max(pv(:,1)), max(pv(:,2))]; % Beschraenke das Gebiet
[vert,tri] = distmesh2d(@dpoly,@huniform,0.1,bbox,false,pv,pv); % Gitter erzeugen

dists1 = zeros(3*length(tri),1); % Vektor mit allen Kantenlaengen initialisieren
miin1 = zeros(length(tri),1); % Vektor fuer alle minimalen Seitenflaechen
maax1 = zeros(length(tri),1); % Vektor fuer alle maximalen Seitenflaechen
counter = 1;
counteer = 1;
for i = 1:length(tri) % Iteriere ueber alle Elemente
    x = vert(tri(i,:),1); % x-Koordinaten des Elements
    y = vert(tri(i,:),2); % y-Koordinaten des Elements
    dists1(counter:counter+2) = sqrt(abs((x-circshift(x,1)).^2+(y-circshift(y,1)).^2));
    % Bestimme die Laenge der Kanten ueber die euklidische Norm
    maax1(counteer) = max(dists1(counter:counter+2));
    miin1(counteer) = min(dists1(counter:counter+2));
    counter = counter + 3;
    counteer = counteer +1;
end
difSeitenlaenge1 = maax1./miin1; % Bestimme das maximale Seitenverhaeltnis 
% fuer jedes Element der Triangulierung
max1 = max(difSeitenlaenge1); % Bestimme das maximale Seitenverhaeltnis 
% ueber alle Elemente
min1 = min(difSeitenlaenge1); % Bestimme das minimale Seitenverhaeltnis 
% ueber alle Elemente

%% Polygon 2
pv = [ 0,  0;...
       0, 22;...
      24, 30;...
      24, 22;...
       0,  0]; % Eckpunkte definieren
bbox = [min(pv(:,1)), min(pv(:,2)); max(pv(:,1)), max(pv(:,2))]; % Beschraenke das Gebiet
[vert2,tri2] = distmesh2d(@dpoly,@huniform,9,bbox,true,pv,pv); % Gitter erzeugen

dists2 = zeros(3*length(tri2),1); % Vektor mit allen Kantenlaengen initialisieren
miin2 = zeros(length(tri),1); % Vektor fuer alle minimalen Seitenflaechen
maax2 = zeros(length(tri),1); % Vektor fuer alle maximalen Seitenflaechen
counter = 1;
counteer = 1;
for i = 1:length(tri2) % Iteriere ueber alle Elemente
    x = vert(tri2(i,:),1); % x-Koordinaten des Elements
    y = vert(tri2(i,:),2); % y-Koordinaten des Elements
    dists2(counter:counter+2) = sqrt(abs((x-circshift(x,1)).^2+(y-circshift(y,1)).^2));
    % Bestimme die Laenge der Kanten ueber die euklidische Norm
    maax2(counteer) = max(dists1(counter:counter+2));
    miin2(counteer) = min(dists1(counter:counter+2));
    counter = counter + 3;
    counteer = counteer +1;
end
difSeitenlaenge2 = maax2./miin2; % Bestimme das maximale Seitenverhaeltnis 
% fuer jedes Element der Triangulierung
max2 = max(difSeitenlaenge2); % Bestimme das maximale Seitenverhaeltnis 
% ueber alle Elemente
min2 = min(difSeitenlaenge2); % Bestimme das minimale Seitenverhaeltnis 
% ueber alle Elemente


%% Dirichletknoten hinzufuegen und plotten
% Toleranz fuer Dirichletknoten einbauen
dirichlet_tol = 10^(-8);
dirichlet = (vert(:,1) >= -0.4 - dirichlet_tol & vert(:,1) <= -0.4 + dirichlet_tol); % Dirichletrand, logischer Vektor

grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen. 
% Macht die uebergabe einfacher und dient als logische Einheit
plotGridDirichlet(grid,1,[],"Triangulierung der Ordnung 1");

%% Testen
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion
order = 1;    %Grad der Basisfunktionen festlegen

[U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

figSolution = plotVectorfieldSolution(vert,tri,U,V,1);

%% Deformierte Flaeche darstellen
[~,~,defVert] = plotDeformationVectors(vert,U,V,order,1);
plotDeformationPolygons(vert,tri,defVert,order,1);
