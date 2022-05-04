function [fig,t] = plotVectorfieldSolution(vert,tri,U,V,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Loesung des Elastizitaetproblems",'NumberTitle','off');
    t = tiledlayout(2,nComparisons,'TileSpacing','Compact','Padding','Compact',TileIndexing='rowmajor');
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end
if isempty(vert)
    return
end
nexttile
trisurf(tri(:,1:3),vert(:,1),vert(:,2),U,'EdgeColor','none');
xlabel('x_1',FontWeight='bold');   ylabel('x_2',FontWeight='bold');
zlabel('\Delta u_1',FontWeight='bold')
title("Deformation in x_1 Richtung")

nexttile
trisurf(tri(:,1:3),vert(:,1),vert(:,2),V,'EdgeColor','none');
xlabel('x_1',FontWeight='bold');    ylabel('x_2',FontWeight='bold');
zlabel('\Delta u_2',FontWeight='bold')
title("Deformation in x_2 Richtung")
end

