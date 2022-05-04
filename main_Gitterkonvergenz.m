clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
addpath('libs/distmesh') % Meshfunktion laden

%% Teste verschiedene Gitterfeinheiten
hVec = 1./(2.^(4:8)); % 1/16, 1/32, 1/64, 1/128, 1/256
maxOrder = 2; % Teste verschiedene Ordnungen der Elemente

verts=cell(length(hVec),maxOrder); tris=cell(length(hVec),maxOrder);
Ucell=cell(length(hVec),maxOrder); Vcell=cell(length(hVec),maxOrder);
Uref=cell(maxOrder,1); Vref=cell(maxOrder,1);
for j = 1:length(hVec)
    for order = 1:maxOrder
        % Gittergenerierung
        [verts{j,order},tris{j,order}] = genMeshSquare(1,1/hVec(j));  % Triangulierung mit Eckknoten erstellen
        [verts{j,order},tris{j,order}] = extendGridLagr(verts{j,order},tris{j,order},order); % Fuer hoehere Ordnung als P1: Hinzufuegen von Knoten
        dirichlet = (verts{j,order}(:,1) == 0); % Dirichletrand, logischer Vektor
        grid = struct("vert",verts{j,order},"tri",tris{j,order},"dirichlet",dirichlet); % Gitter in eine Structure  bringen
        
        % PDE
        E = 210; nu = 0.3; % Materialparameter
        f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
        gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
        
        % Aufstellen der Loesung
        [Ucell{j,order},Vcell{j,order}] = elastSolver(grid,E,nu,f,gD,order);
        
        % Die Loesungen der verschiedenen Gitter anpassen, sodass sie nur die
        % Werte an den uebereinstimmenden Punkten enthalten, also die der groebsten
        % Triangulierung
        % sonst: spaeter Dimensionsprobleme beim Vergleich mit der Referenzloesung
        if (j > 1)
            % Pruefe, ob Knoten in  beiden Knotenlisten enthalten sind
            % ind: logischer Vektor; loc: gibt Position des Knotens in der anderen Liste an
            [ind,loc] = ismember(verts{j,order},verts{1,order},'rows'); 
            
            % Eliminiere Werte an nicht doppelt vorkommenden Knoten
            Ucell{j,order}=Ucell{j,order}(ind,:); 
            Vcell{j,order}=Vcell{j,order}(ind,:);
            
            % Sortiere Werte um, sodass die Knotennummerierungen uebereinstimmen (notwendig fuer Elemente hoeherer Ordnung)
            Ucell{j,order}(loc(ind),:)=Ucell{j,order}; 
            Vcell{j,order}(loc(ind),:)=Vcell{j,order};
        end
        
        % Die Loesung auf dem feinsten Gitter dient als Referenzloesung
        if (j == length(hVec))
            Uref{order}=Ucell{end,order};
            Vref{order}=Vcell{end,order};
        end
    end
end

% Berechnung der Abweichung der Loesungen von der Referenzloesung
Udiff=cell(length(hVec)-1,maxOrder); Vdiff=cell(length(hVec)-1,maxOrder);
for j=1:length(hVec)-1
    for order = 1:maxOrder
        % Messung der Abweichung ueber gemittelte Zeilensummennorm
        Udiff{j,order}=norm((Ucell{j,order}-Uref{order})',1)/length(Ucell{j,order});
        Vdiff{j,order}=norm((Vcell{j,order}-Vref{order})',1)/length(Vcell{j,order});
    end
end
Udiff=cell2mat(Udiff); Vdiff=cell2mat(Vdiff); % Cell-Array in Matrix umwandeln

% Plotten der Ergebnisse
figure("Name","Gitterkonvergenz")
for order = 1:maxOrder
    subplot(1,maxOrder,order)
    plot(hVec(1:end-1),Udiff(:,order),'blue'); hold on;
    plot(hVec(1:end-1),Vdiff(:,order),'red'); hold on;
    scatter(hVec(1:end-1),Udiff(:,order),'filled','blue'); hold on;
    scatter(hVec(1:end-1),Vdiff(:,order),'filled','red');
    xlim([0,hVec(1)])
    ylim([0,max(max(Udiff(:,order)),max(Vdiff(:,order)))])
    xlabel('Gitterfeinheit',FontWeight='bold')
    ylabel('Abweichung in gemittelter Zeilensummennorm',FontWeight='bold')
    legend('Abweichung in x1-Richtung','Abweichung in x2-Richtung','Location','northwest')
    title('Abweichung der Loesung von der Referenzloesung')
    subtitle(sprintf('Ordnung %g',order))
end

