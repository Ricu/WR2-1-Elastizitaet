function [fig,t] = plotVectorfieldSolution(vert,tri,U,V,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Loesung des Elastizitaetproblems",'NumberTitle','off');
    t = tiledlayout(1,nComparisons,TileIndexing='rowmajor');
else
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

nexttile
trisurf(tri(:,1:3),vert(:,1),vert(:,2),U);
title("U: x_1 Richtung")
nexttile
trisurf(tri(:,1:3),vert(:,1),vert(:,2),V);
title("V: x_2 Richtung")
end

