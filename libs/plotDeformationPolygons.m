function [fig,t] = plotDeformationPolygons(tri,vert,defVert,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformierung",'NumberTitle','off');
    t = tiledlayout(2,nComparisons);
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

nexttile
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
axis equal tight;
nexttile
patch('vertices',defVert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
axis equal tight;
end

