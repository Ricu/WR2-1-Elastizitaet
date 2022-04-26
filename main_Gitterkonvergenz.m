clear; clc; % Konsolen Output und Variablen loeschen
addpath('libs') % Hilfsfunktionen laden
%% Teste verschiedene Gitterfeinheiten
maxOrder = 2;
hVec = 1./(2.^(4:7)); % 16, 32, 64, 128
nu = 0.3;

Ucell=cell(length(hVec),maxOrder); Vcell=cell(length(hVec),maxOrder);

for j = 1:length(hVec)
    h = hVec(j);
    for order = 1:maxOrder
        [vert,tri] = genMeshSquare(1,1/h);
        [vert,tri] = extendGridLagr(vert,tri,order);
        dirichlet = (vert(:,1) == 0); % Dirichletrand, logischer Vektor
        grid = struct("vert",vert,"tri",tri,"dirichlet",dirichlet); % Gitter in eine Structure  bringen.
        
        % PDE
        E = 210;  % Materialparameter
        f = @(x,y) [ones(size(x));ones(size(y))]; % Volumenkraft
        gD = @(x) 0*x; % Dirichlet-Randwertfunktion, x=[x_1;x_2]
        
        [Ucell{j,order},Vcell{j,order}] = elastSolver(grid,E,nu,f,gD,order); % Problem loesen      
    end
end

Uref=cell(maxOrder,1); Vref=cell(maxOrder,1);
for order=1:maxOrder
    Uref{order}=Ucell{end,order};
    Vref{order}=Vcell{end,order};
end

Udiff=cell(length(hVec)-1,maxOrder); Vdiff=cell(length(hVec)-1,maxOrder);
for j=1:length(hVec)-1
    for order = 1:maxOrder
        Udiff{j,order}=abs(norm(Ucell{j,order},1)/size(Ucell{j,order},1)-norm(Uref{order},1)/size(Uref{order},1));
        Vdiff{j,order}=abs(norm(Vcell{j,order},1)/size(Vcell{j,order},1)-norm(Vref{order},1)/size(Vref{order},1));
    end
end
Udiff=cell2mat(Udiff); Vdiff=cell2mat(Vdiff);

figure("Name","Gitterkonvergenz")
for order = 1:maxOrder
    subplot(1,maxOrder,order)
    plot(hVec(1:end-1),Udiff(:,order));
    hold on;
    plot(hVec(1:end-1),Vdiff(:,order));
    xlabel('Gitterfeinheit')
    legend('Abweichung in x1-Richtung','Abweichung in x2-Richtung')
    title('Abweichung der Loesung von der Referenzloesung fuer verschiedene Gitterfeinheiten')
    subtitle(sprintf('Ordnung %g',order))
end

