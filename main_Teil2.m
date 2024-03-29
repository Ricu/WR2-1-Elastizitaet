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
h0 = 0.1; % Angestrebte Kantenlaenge
[vert,tri] = distmesh2d(@dpoly,@huniform,h0,bbox,true,pv,1,pv); % Gitter erzeugen mit Gittergenerator

edgeLengthsPol1 = zeros(3*length(tri),1); % Vektor mit allen Kantenlaengen initialisieren
counter = 1;
for i = 1:length(tri) % Iteriere ueber alle Elemente
    x = vert(tri(i,:),1); % x-Koordinaten des Elements
    y = vert(tri(i,:),2); % y-Koordinaten des Elements
    % Bestimme die Kantenlaenge mit euklidischer Norm
    edgeLengthsPol1(counter:counter+2) = sqrt(abs((x-circshift(x,1)).^2+(y-circshift(y,1)).^2));
    counter = counter + 3;
end
% Entferne doppelt vorkommende Kanten aus der Liste
edgeLengthsPol1 = dropDuplicateEdges(edgeLengthsPol1,tri);

% Kantenlaengenverteilung plotten
figure('Name',"Kantenlaengenverteilung fuer Polygon 1",NumberTitle="off")
histogram(edgeLengthsPol1)
xline(h0, '-r','h_0','LineWidth',2);
ylabel("Frequenz",FontWeight='bold'); xlabel("Kantenlaenge",FontWeight='bold');
title(sprintf("Kantenlaengenverteilung fuer Polygon 1 mit h_0 = %.1f", h0))

%% Polygon 2
pv = [ 0,  0;...
       0, 22;...
      24, 30;...
      24, 22;...
       0,  0]; % Eckpunkte definieren
bbox = [min(pv(:,1)), min(pv(:,2)); max(pv(:,1)), max(pv(:,2))]; % Beschraenke das Gebiet
h0 = 6; % Angestrebte Kantenlaenge
[vert2,tri2] = distmesh2d(@dpoly,@huniform,h0,bbox,true,pv,2,pv); % Gitter erzeugen mit Gittergenerator

edgeLengthsPol2 = zeros(3*length(tri2),1); % Vektor mit allen Kantenlaengen initialisieren
counter = 1;
for i = 1:length(tri2) % Iteriere ueber alle Elemente
    x = vert2(tri2(i,:),1); % x-Koordinaten des Elements
    y = vert2(tri2(i,:),2); % y-Koordinaten des Elements
    % Bestimme die Laenge der Kanten ueber die euklidische Norm
    edgeLengthsPol2(counter:counter+2) = sqrt(abs((x-circshift(x,1)).^2+(y-circshift(y,1)).^2));
    counter = counter + 3;
end
% Entferne doppelt vorkommende Kanten aus der Liste
edgeLengthsPol2 = dropDuplicateEdges(edgeLengthsPol2,tri2);

% Kantenlaengenverteilung plotten
figure('Name',"Kantenlaengenverteilung fuer Polygon 2",NumberTitle="off")
histogram(edgeLengthsPol2, [0, 0.7*h0, 0.9*h0,1.1*h0, 1.3*h0, 2*h0])
xline(h0, '-r','h_0','LineWidth',2);
ylabel("Frequenz",FontWeight='bold'); xlabel("Kantenlaenge",FontWeight='bold');
title(sprintf("Kantenlaengenverteilung fuer Polygon 2 mit h_0 = %.1f", h0))

%% Betrachte ab hier nur Polygon 1
%% Dirichletknoten hinzufuegen und Triangulierung plotten
% Toleranz fuer Dirichletknoten einbauen
dirichlet_tol = 10^(-8);
% Dirichletrand, logischer Vektor
xdirichlet=-0.4;
dirichlet = (vert(:,1) >= xdirichlet - dirichlet_tol & vert(:,1) <= xdirichlet + dirichlet_tol); 
grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Struktur bringen 
plotGridDirichlet(grid,1,[],"Triangulierung der Ordnung 1");

%% PDE
E = 210; nu = 0.3; % Materialparameter
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion
order = 1;  % Grad der Basisfunktionen festlegen

%% Problem loesen
[U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

%% Loesung plotten
figSolution = plotVectorfieldSolution(vert,tri,U,V,1);

%% Deformierte Flaeche darstellen
% Plot der veraenderten Position der Knoten mit Verschiebungsvektoren
[~,~,defVert] = plotDeformationVectors(vert,U,V,order,1);
% Plot des defomierten Polygons
plotDeformationPolygons(vert,tri,defVert,order,1);

%% Hilfsfunktion
% Entfernt doppelt vorkommende Kanten aus der Liste
function edgeLengthList_reduced = dropDuplicateEdges(edgeLengthList,tri)
tri_circ = circshift(tri,1,2);
edges = [tri(:), tri_circ(:)];
edges = sort(edges,2);
[~, ia, ~] = unique(edges,'rows');
edgeLengthList_reduced = edgeLengthList(ia);
end