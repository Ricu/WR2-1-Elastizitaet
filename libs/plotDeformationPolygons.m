function [fig,t] = plotDeformationPolygons(vert,tri,defVert,order,nComparisons,fig,customTitle)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
    t = tiledlayout(1,nComparisons,'TileSpacing','Compact','Padding','Compact',TileIndexing='rowmajor');
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
    xlabel(t,'x_1',FontWeight='bold'); ylabel(t,'x_2',FontWeight='bold');
end
if isempty(vert)
    return
end
if ~(exist('customTitle','var')) || isempty(customTitle)
    customTitle = "Polygon nach Deformation";
end
% Knotennummerierung fuer patch aendern; Nummerierung gegen den
% Uhrzeigersinn
if order == 2
    tri = tri(:,[1,5,2,6,3,4]);
end

nexttile
patch('vertices',defVert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
title(customTitle)
axis equal tight;

end

