function [fig,t] = plotDeformationPolygons(vert,tri,defVert,order,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
    t = tiledlayout(1,nComparisons,TileIndexing='rowmajor');
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end
if isempty(vert)
    return
end
% Knotennummerierung fuer patch aendern; Nummerierung gegen den
% Uhrzeigersinn
if order == 2
    tri = tri(:,[1,5,2,6,3,4]);
end

nexttile
patch('vertices',defVert,'faces',tri,'edgecol','k','facecol',[.8,.9,1]);
subtitle("Polygon nach Deformation")
axis equal tight;
xlabel('x_1',FontWeight='bold'); ylabel('x_2',FontWeight='bold');
end

