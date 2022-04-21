function [fig,t] = plotDeformationVectors(order,vert,defVert,U,V,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
    t = tiledlayout(1,nComparisons);
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end
xmin = min(min(vert(:,1)),min(defVert(:,1)));
xmax = max(max(vert(:,1)),max(defVert(:,1)));
ymin = min(min(vert(:,2)),min(defVert(:,2)));
ymax = max(max(vert(:,2)),max(defVert(:,2)));
xdist = xmax - xmin;
ydist = ymax - ymin;
xmin = xmin - 0.05 * xdist; xmax = xmax + 0.05 * xdist;
ymin = ymin - 0.05 * ydist; ymax = ymax + 0.05 * ydist;
nexttile
scatter(vert(:,1),vert(:,2),'k','filled'); hold on; % Urspruengliche Flaeche plotten
scatter(defVert(:,1),defVert(:,2),46); % Deformierte Flaeche plotten
quiver(vert(:,1),vert(:,2),U,V,0) % Berechnetes Vektorfeld (Verschiebung) plotten
legend("Original","Deformiert","Verschiebung") % Legende hinzufuegen
title(sprintf("Ordnung %g",order))
xlim([xmin,xmax])
ylim([ymin,ymax])
end

