function [fig,t] = plotDeformationPolygons(tri,vert,order,defVert,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
    t = tiledlayout(2,nComparisons,TileIndexing='rowmajor');
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

% Knotennummerierung fuer patch aendern; Nummerierung gegen den
% Uhrzeigersinn
if order == 2
    tri_copy=tri;
    tri(:,1)=tri_copy(:,1);
    tri(:,2)=tri_copy(:,5);
    tri(:,3)=tri_copy(:,2);
    tri(:,4)=tri_copy(:,6);
    tri(:,5)=tri_copy(:,3);
    tri(:,6)=tri_copy(:,4);
end

nexttile
patch('vertices',vert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
title(sprintf("Ordnung %g",order))
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subtitle("Polygon vor Deformation")
axis equal tight;
nexttile
patch('vertices',defVert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
subtitle("Polygon nach Deformation")
axis equal tight;
end

