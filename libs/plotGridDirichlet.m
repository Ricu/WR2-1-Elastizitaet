function [fig,t] = plotGridDirichlet(grid,nComparisons,fig,customTitle)

if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Triangulierung",'NumberTitle','off');
    t = tiledlayout(1,nComparisons);
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

if ~(exist('customTitle','var')) || isempty(customTitle)
    customTitle = "Triangulierung entsprechender Ordnung";
end

if isempty(grid)
    return
end
vert = grid.vert; tri = grid.tri; dirichlet = grid.dirichlet; % Gitterdaten extrahieren

nexttile
patch('vertices',vert,'faces',tri(:,1:3),'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten
hold on; 
axis equal tight
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
title(customTitle)
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufuegen
end

