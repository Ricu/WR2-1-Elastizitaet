clear; clc;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

% In diesem Skript vergleichen wir die von dem PCG Verfahren benoetigte
% Anzahl an Iterationen zum Loesen des Galerkin System fuer lineare
% Elastizitaet. Es werden P1 und P2 Elemente verglichen.

%% Verschiedene Werte fuer die Querkontraktionszahl (nu->1/2) bei verschiedenen Ordnungen testen
maxOrder = 2;
nuVec = [ 0.45, 0.49, 0.4999 , 0.499999]; % [0.05, 0.2, 0.4, 0.45, 0.49];

%% PDE
E = 210; nu = 0.3; % Materialparameter
% E = 21*10^5; nu = 0.28; % Stahl
% E = 0.037; nu = 0.485; % Radiergummi
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]

%% Figure
figDefPol = plotDeformationPolygons([],[],[],[],length(nuVec)*maxOrder);

numIter=zeros(length(nuVec),maxOrder);
for order=1:maxOrder
    for j=1:length(nuVec)
        %% Mesh
        h = 1/16; % Gitterfeinheit
        [vert,tri] = genMeshSquare(1,1/h); % Punkteliste und Elementeliste erstellen
        [vert,tri] = extendGridLagr(vert,tri,order); % Fuer hoehere Ordnung als P1: Hinzufuegen von Knoten
        dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
        dirichlet2 = repelem(dirichlet,2);  % Dirichletrand um die 2. Komponente erweitern
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen.
        
        %% Loese das System mit iterativem Loeser
        [~,~,K,F] = elastSolver(grid,E,nuVec(j),f,gD,order); % Steifigkeitsmatrix und Lastvektor aufstellen
        VK = eye(size(K)); % Definiere einen Vorkonditionierer
        [x,numIter(j,order)] = PCG(K,F,VK); % Loese das Sytem mit PCG
        
        d = zeros(length(F),1); % Loesungsvektor initialisieren
        d(dirichlet2) = gD(vert(dirichlet2)); % Dirichletwerte eintragen
        d(~dirichlet2) = x; % Systemloesung an richtige Stellen schreiben
        U = d(1:2:end); % Loesung in x_1 Richtung extrahieren
        V = d(2:2:end); % Loesung in x_2 Richtung extrahieren
        
         %% Deformierte Flaeche darstellen
        defVert = vert; % Deformierte Liste initialisieren
        defVert(:,1) = defVert(:,1) + U; % Deformierung in x_1 Richtung
        defVert(:,2) = defVert(:,2) + V; % Deformierung in x_2 Richtung
        % Plot des defomierten Polygons
        figDefPol = plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol,sprintf("Ordnung=%g, nu=%g",order,nuVec(j)));
    end
end

%% Konvergenzvergleich plotten
figure("Name","Iterativer Loeser")
plot(nuVec,numIter(:,1),'blue'); hold on;
plot(nuVec,numIter(:,2),'red'); hold on;
scatter(nuVec,numIter(:,1),'filled','blue'); hold on;
scatter(nuVec,numIter(:,2),'filled','red');
xlabel('Materialparameter nu',FontWeight='bold')
ylabel('Anzahl Iterationen',FontWeight='bold')
legend('P1 Elemente','P2 Elemente','Location','northwest')
title('Anzahl Iterationen von PCG zur Loesung des Systems')



