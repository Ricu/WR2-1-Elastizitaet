clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden
%% Teste verschiedene Materialparameter
maxOrder = 1;
h = 1/16;
nuVec = [0.1, 0.3, 0.49]; % [0.05, 0.2, 0.4, 0.45, 0.49];
EVec = [100, 210, 500];
order = 1;

pv = [ 0  , 0  ;...
       4  , 0  ;...
       4  , 1  ;...
       0  , 1  ;...
       0  , 0  ]; % Eckpunkte definieren

bbox = [min(pv(:,1)), min(pv(:,2)); max(pv(:,1)), max(pv(:,2))];
[vert,tri] = distmesh2d(@dpoly,@huniform,0.4,bbox,true,pv,pv);
dirichlet_tol = 10^(-8);
xdirichlet = 0;
dirichlet = (vert(:,1) >= xdirichlet - dirichlet_tol & vert(:,1) <= xdirichlet + dirichlet_tol); % Dirichletrand, logischer Vektor

figDefPol = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
t = tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
xlabel(t,'x_1',FontWeight='bold')
ylabel(t,'x_2',FontWeight='bold')
nexttile(1,[1 3])
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[1,1,1]);
ylim([-0.5,1.5]); xlim([0,4.5]);
axes = cell(length(EVec),length(nuVec));

for i = 1:length(EVec)
    E = EVec(i);
    for j = 1:length(nuVec)
        nu = nuVec(j);
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen.

        % PDE
        f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
        gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]

        [U,V] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen

        deformed_area = vert; % Deformierte Liste initialisieren
        deformed_area(:,1) = deformed_area(:,1) + U; % Deformierung in x_1 Richtung
        deformed_area(:,2) = deformed_area(:,2) + V; % Deformierung in x_2 Richtung

        %% Plot        
        figure(figDefPol)
        t = get(figDefPol,'children');
        axes{i,j} = nexttile;
        patch('vertices',deformed_area,'faces',tri,'edgecol','k','facecol',[1,1,1]);
        title(sprintf("nu = %3.2f, E = %i",nu,E))
        axis equal tight; 
   end
end
linkaxes([axes{:}])
