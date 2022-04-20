function [fig,t] = plotDeformationVectors(order,vert,defVert,U,V,nComparisons,fig)
if ~(exist('fig','var')) || isempty(fig)
    fig = figure("Name","Gebietsvergleich: vor und nach Deformation",'NumberTitle','off');
    t = tiledlayout(1,nComparisons);
else
    figure(fig)
    t = get(fig,'children'); % Extrahiere das TiledLayoutObject unter der Annahme, das alle Plots mit diesem erstellt wurden
end

nexttile
scatter(vert(:,1),vert(:,2),'k','filled'); hold on; % Urspruengliche Flaeche plotten
scatter(defVert(:,1),defVert(:,2),46); % Deformierte Flaeche plotten
quiver(vert(:,1),vert(:,2),U,V,0) % Berechnetes Vektorfeld (Verschiebung) plotten
legend("Original","Deformiert","Verschiebung") % Legende hinzufuegen
title(sprintf("Ordnung %g",order))
xlim([0,1.1])
ylim([0,1.1])
end

