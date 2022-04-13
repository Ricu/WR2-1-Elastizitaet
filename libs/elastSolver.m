function [U,V] = elastSolver(grid,E,nu,force,gD)

vert = grid.vert; tri = grid.tri; dirichlet = grid.dirichlet; % Gitterdaten extrahieren
% Da jeder Knoten 2 Basisfunktionen verwendet, muss dies auch im logischen
% Vektor berücksichtigt werden
ind = [2*find(dirichlet)-1 ; 2*find(dirichlet)]; % Dirichletrand um die 2. Komponente erweitern
dirichlet2 = false(2*length(vert),1); % Neuen logischen Vektor erstellen
dirichlet2(ind) = true; % Dirichletknoten markieren

[mu,lambda]=enu2lame(E,nu); % Materialparameter extrahieren
[K,~,F] = elastAssemble(vert',tri',mu,lambda,force,dirichlet2,gD); % Assemblierungsroutine aufrufen


d = zeros(length(F),1); % Lösungsvektor initialisieren
d(dirichlet2) = gD(vert(dirichlet2)); % Dirichletwerte eintrage
d(~dirichlet2) = K\F; % Galerkin System lösen
U = d(1:2:end); % Lösung in x_1 Richtung extrahieren
V = d(2:2:end); % Lösung in x_2 Richtung extrahieren
end

function [K,M,F,K_full,F_full] = elastAssemble(p,t,lambda,mu,force,dirichlet,gD)
% Input: p als matrix mit allen Knoten
% Input: e als matrix mit allen Kanten
% Input: t als matrix mit allen Verbindungen der Trainagulierung

% Output: K als global assemblierte Steifigkeitsmatrix
% Output: M als global assemblierte Massenmatrix
% Output: F als global assemblierter load vector
ndof=2*size(p,2); % absolute Anzahl der Freiheitsgrade entspricht ...
% der Anzahl an Knoten*2
% K=sparse(ndof,ndof); % initialisiere Steifigkeitsmatrix
% M=sparse(ndof,ndof); % initialisiere Massenmatrix
F=zeros(ndof,1); % initialisiere Ladungsvektor

% Matrizen vorbereiten, in die die Index-Wertepaare für K&M gespeichert
% werden
nBaseFun = 6; % Anzahl Basisfunktion für Order = 1: 6
% Anzahl Basisfunktion für Order (Lagrange) = k: (k+2)*(k+1);
nEleLoc = nBaseFun^2; % Anzahl Elemente der lokalen Steifigkeitsmatrix
K_val = zeros(nEleLoc*length(t),1);
M_val = K_val;
iIndex = K_val;
jIndex = K_val;

dofs=zeros(6,1); % Initialisiere Vektor, in dem die Knotenindizes gespeichert werden
for i=1:size(t,2) % Über die Elemente iterieren
    nodes=t(1:3,i); % Indizes der zum aktuellen Element gehörenden physikalischen Knoten Indizes
    x=p(1,nodes); y=p(2,nodes); % Koordinaten der Knoten
    
    f=force(x,y); % Volumenkraft an aktuellen Knoten auswerten
    KK=elasticStiffness(x,y,lambda,mu); % lokale Steifigkeitsmatrix berechnen
    MK=elasticMass(x,y); % lokale Massenmatrix berechnen
    fK=[f(1,1) f(2,1) f(1,2) f(2,2) f(1,3) f(2,3)]'; % Volumenkraft an den aktuellen Knoten bestimmen
    FK=MK*fK; % lokalen Ladungsvektor bestimmen
    
    % Assembliere durch addieren der lokalen Matrizen an die richtigen
    % Stellen der globalen Matrizen
    dofs(2:2:end)=2*nodes;          % Knoten Indizes um die 2. Komponente erweitern
    dofs(1:2:end)=2*nodes-1;        % Knoten Indizes um die 2. Komponente erweitern
    %     K(dofs,dofs)=K(dofs,dofs)+KK;   % addieren zur Steifigkietsmatrix
    %     M(dofs,dofs)=M(dofs,dofs)+MK;   % addieren zur Massenmatrix
    %     F(dofs)=F(dofs)+FK;             % addieren zum load vector
    
    
    iIndex((i-1)*nEleLoc+1:i*nEleLoc) = reshape(repmat(dofs',nBaseFun,1),nEleLoc,1); % i-Indizes speichern
    jIndex((i-1)*nEleLoc+1:i*nEleLoc) = repmat(dofs,nBaseFun,1); % j-Indizes speichern
    K_val((i-1)*nEleLoc+1:i*nEleLoc) = reshape(KK,nEleLoc,1); % Werte für die Steifigkeitsmatrix speichern
    M_val((i-1)*nEleLoc+1:i*nEleLoc) = reshape(MK,nEleLoc,1); % Werte für die Massenmatrix speichern
    F(dofs) = F(dofs) + FK; % Der Lastvektor kann ohne Umwege aktualisiert werden
end
K = sparse(iIndex,jIndex,K_val,ndof,ndof); % Anhand der zuvor erstellten Listen die Steifigkeitsmatrix erstellen
M = sparse(iIndex,jIndex,M_val,ndof,ndof); % Anhand der zuvor erstellten Listen die Massenmatrix erstellen


K_full = K; % Speichere die Steifigkeitsmatrix inklusive der Dirichletknoten
F_full = F - K(:,dirichlet)* gD(p(dirichlet)); % Speichere den Lastvektor inklusive der Dirichletknoten

K = K(~dirichlet,~dirichlet); % Dirichletknoten aus der Steifigkeitsmatrix eliminieren
F = F(~dirichlet)- K(:,dirichlet)*gD(p(dirichlet)); % Dirichletknoten aus dem Lastvektor eliminieren
% TODO VAL, muss G_D hier abgezogen werden? Glaube ja.
end


function KK = elasticStiffness(x,y,mu,lambda)
% Input: x und y Koordinaten eines Elements
% Input: Lame Parameter mu und lambda

% Output: Lokale Steifigkeitsmatrix

[area,b,c]=hatGradients(x,y); % Fläche des Elements und Gradienten der Basisfkt. bestimmen
D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0]; % Elastizität Matrix aufstellen
% Spannungsmatrix aufstellen
BK=[b(1) 0 b(2) 0 b(3) 0 ;
    0 c(1) 0 c(2) 0 c(3);
    c(1) b(1) c(2) b(2) c(3) b(3)];
KK=BK'*D*BK*area; % Lokale Steifigkeitsmatrix bestimmen
end

function MK = elasticMass(x,y)
area=polyarea(x,y); % Fläche mittels polyarea bestimmen
% lokale Massenmatrix bestimmen
MK=[2 0 1 0 1 0;
    0 2 0 1 0 1;
    1 0 2 0 1 0;
    0 1 0 2 0 1;
    1 0 1 0 2 0;
    0 1 0 1 0 2]*area/12;
end

function [mu,lambda] = enu2lame(E,nu)
% Input: E ist der Young modulus
% Input: nu ist die Poisson ratio

% Output: Lame Parameter mu und lambda
mu=E/(2*(1+nu)); % Mu berechnen
lambda=E*nu/((1+nu)*(1-2*nu)); % Lambda berechnen
end

function [area,b,c] = hatGradients(x,y)
% Input: x und y Koordinaten eines Elements

% Output: ...
area=polyarea(x,y); % Fläche mittels polyarea bestimmen
b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area; % Gradient bestimmen
c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area; % Gradient bestimmen
end