                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

%% Teste verschiedene Werte für nu und verschiedene Ordnungen
maxOrder = 2;
nuVec = [0.3, 0.45, 0.49]; % [0.05, 0.2, 0.4, 0.45, 0.49];

%% PDE
E = 210; nu = 0.3; % Materialparameter
% E = 21*10^5; nu = 0.28; % Stahl
% E = 0.037; nu = 0.485; % Radiergummi
f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]

%% Figure
figDefPol = plotDeformationPolygons([],[],[],[],length(nuVec)*maxOrder);

for order=1:maxOrder
    for j=1:length(nuVec)
        %% Mesh
        h = 1/16;
        [vert,tri] = genMeshSquare(1,1/h); % Punkteliste und Elementeliste erstellen
        [vert,tri] = extendGridLagr(vert,tri,order); % Gittererweiterung bei hoeherer Ordnung
        dirichlet = (vert(:,1) == 0 | vert(:,2) == 0 | vert(:,1) == 1 | vert(:,2) == 1); % Dirichletrand, logischer Vektor
        dirichlet2 = repelem(dirichlet,2);  % Erweitere den Dirichletrand um die 2. Komponente
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen.
        
        %% Loese das System mit iterativem Loeser
        [~,~,K,F] = elastSolver(grid,E,nuVec(j),f,gD,order); % Steifigkeitsmatrix und Lastvektor aufstellen
        VK = eye(size(K)); % Definiere einen Vorkonditionierer
        [x,k] = PCG(K,F,VK); % Loese das Sytem mit PCG
        
        d = zeros(length(F),1); % Loesungsvektor initialisieren
        d(dirichlet2) = gD(vert(dirichlet2)); % Dirichletwerte eintragen
        d(~dirichlet2) = x; % Systemloesung an richtige Stellen schreiben
        U = d(1:2:end); % Loesung in x_1 Richtung extrahieren
        V = d(2:2:end); % Loesung in x_2 Richtung extrahieren
        
        defVert = vert; % Deformierte Liste initialisieren
        defVert(:,1) = defVert(:,1) + U; % Deformierung in x_1 Richtung
        defVert(:,2) = defVert(:,2) + V; % Deformierung in x_2 Richtung
        
        % Deformierte Flaeche darstellen
        figDefPol = plotDeformationPolygons(vert,tri,defVert,order,2,figDefPol);
    end
end
