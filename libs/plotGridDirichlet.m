function [fig,t] = plotGridDirichlet(vert,tri,dirichlet,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure();
    t = tiledlayout(1,nComparisons);
else
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

nexttile
patch('vertices',vert,'faces',tri(:,1:3),'edgecol','k','facecol',[.8,.9,1]); % Triangulierung plotten %TODO VALE triplot nutzen
hold on; 
scatter(vert(dirichlet,1),vert(dirichlet,2),[],"r") % Dirichletknoten markieren
legend("Triangulierung","Dirichletrand Knoten") % Legende hinzufuegen
end

