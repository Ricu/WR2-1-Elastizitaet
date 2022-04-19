function [fig,t] = plotDeformationPolygons(vert,defVert,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure();
    t = tiledlayout(1,nComparisons);
else
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

nexttile
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
axis equal tight;
nexttile
patch('vertices',defVert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
axis equal tight;
end

