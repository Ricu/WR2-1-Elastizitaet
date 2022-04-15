function [U,V] = elastSolver(grid,E,nu,f,gD)

vert = grid.vert; tri = grid.tri; dirichlet = grid.dirichlet; % Gitterdaten extrahieren
% Da jeder Knoten 2 Basisfunktionen verwendet, muss dies auch im logischen
% Vektor berücksichtigt werden
ind = [2*find(dirichlet)-1 ; 2*find(dirichlet)]; % Dirichletrand um die 2. Komponente erweitern %TODO zweite komponente
dirichlet2 = false(2*length(vert),1); % Neuen logischen Vektor erstellen
dirichlet2(ind) = true; % Dirichletknoten markieren

[mu,lambda]=enu2lame(E,nu); % Materialparameter extrahieren
[K,~,F] = elastAssemble(vert',tri',mu,lambda,f,dirichlet2,gD); % Assemblierungsroutine aufrufen


d = zeros(length(F),1); % Lösungsvektor initialisieren
d(dirichlet2) = gD(vert(dirichlet2)); % Dirichletwerte eintrage
d(~dirichlet2) = K\F; % Galerkin System lösen
U = d(1:2:end); % Lösung in x_1 Richtung extrahieren
V = d(2:2:end); % Lösung in x_2 Richtung extrahieren
end

function [K,M,F,K_dir,F_dir] = elastAssemble(p,t,lambda,mu,f,dirichlet,gD)
% Input: p als matrix mit allen Knoten
% Input: e als matrix mit allen Kanten
% Input: t als matrix mit allen Verbindungen der Trainagulierung

% Output: K als global assemblierte Steifigkeitsmatrix
% Output: M als global assemblierte Massenmatrix
% Output: F als global assemblierter load vector

% Load basic functions and quadrature data
[~,d_phihat] = baseFun(1); %TODO: variable order, statt hard coding 1
% if order == 1
%     quad_low = load("quad_formeln.mat").quadratur_P1;
%     quad_high = load("quad_formeln.mat").quadratur_P5;
% elseif order == 2
%     quad_low = load("quad_formeln.mat").quadratur_P2;
%     quad_high = load("quad_formeln.mat").quadratur_P5;
% end

n_nodes = 2*size(p,2); % absolute Anzahl der Freiheitsgrade entspricht ...
% der Anzahl an Knoten*2
F=zeros(n_nodes,1); % initialisiere Ladungsvektor

% Matrizen vorbereiten, in die die Index-Wertepaare für K&M gespeichert
% werden
nBaseFun = 6; % Anzahl Basisfunktion für Order = 1: 6
% Anzahl Basisfunktion für Order (Lagrange) = k: (k+2)*(k+1);
nEleLoc = nBaseFun^2; % Anzahl Elemente der lokalen Steifigkeitsmatrix
K_val = zeros(nEleLoc*length(t),1);
M_val = K_val;
iIndex = K_val;
jIndex = K_val;

node_ind2=zeros(6,1); % Initialisiere Vektor, in dem die Knotenindizes gespeichert werden
for i=1:size(t,2) % Über die Elemente iterieren
    node_ind=t(1:3,i); % die zum aktuellen Element gehörenden physikalischen Knoten Indizes
    x=p(1,node_ind); y=p(2,node_ind); % Koordinaten der Knoten
    
    f_eval=f(x,y); % Volumenkraft an aktuellen Knoten auswerten
    KK=elasticStiffness(x,y,lambda,mu,d_phihat); % lokale Steifigkeitsmatrix berechnen
    MK=elasticMass(x,y); % lokale Massenmatrix berechnen
    fK=[f_eval(1,1) f_eval(2,1) f_eval(1,2) f_eval(2,2) f_eval(1,3) f_eval(2,3)]'; % Volumenkraft an den aktuellen Knoten bestimmen
    FK=MK*fK; % lokalen Lastsvektor bestimmen
    
    % Assembliere durch addieren der lokalen Matrizen an die richtigen
    % Stellen der globalen Matrizen
    node_ind2(2:2:end)=2*node_ind;          % Knoten Indizes um die 2. Komponente erweitern
    node_ind2(1:2:end)=2*node_ind-1;        % Knoten Indizes um die 2. Komponente erweitern
   
    
    iIndex((i-1)*nEleLoc+1:i*nEleLoc) = reshape(repmat(node_ind2',nBaseFun,1),nEleLoc,1); % i-Indizes speichern
    jIndex((i-1)*nEleLoc+1:i*nEleLoc) = repmat(node_ind2,nBaseFun,1); % j-Indizes speichern
    K_val((i-1)*nEleLoc+1:i*nEleLoc) = reshape(KK,nEleLoc,1); % Werte für die Steifigkeitsmatrix speichern
    M_val((i-1)*nEleLoc+1:i*nEleLoc) = reshape(MK,nEleLoc,1); % Werte für die Massenmatrix speichern
    F(node_ind2) = F(node_ind2) + FK; % Der Lastvektor kann ohne Umwege aktualisiert werden
end
K = sparse(iIndex,jIndex,K_val,n_nodes,n_nodes); % Anhand der zuvor erstellten Listen die Steifigkeitsmatrix erstellen
M = sparse(iIndex,jIndex,M_val,n_nodes,n_nodes); % Anhand der zuvor erstellten Listen die Massenmatrix erstellen


K_dir = K; % Speichere die Steifigkeitsmatrix inklusive der Dirichletknoten
F_dir = F - K(:,dirichlet)* gD(p(dirichlet)); % Speichere den Lastvektor inklusive der Dirichletknoten %TODO nicht _full nennen

double_points = zeros(2,2*length(p));
double_points(2:2:end) = p; % TODO Vale: doppelte punkte bessser abspeichern?
double_points(1:2:end) = p;
F = F(~dirichlet)- K(~dirichlet,dirichlet)*gD(double_points(dirichlet)); %
K = K(~dirichlet,~dirichlet); % Dirichletknoten aus der Steifigkeitsmatrix eliminieren
% TODO VAL, muss G_D hier abgezogen werden? Glaube ja.

% Dirichletknoten aus dem Lastvektor eliminieren %TODO
% F = F(~dirichlet)- K(:,dirichlet)*gD(p(dirichlet)); % Dirichletknoten aus dem Lastvektor eliminieren

end


function KK = elasticStiffness(x,y,mu,lambda,d_phihat)
% Input: x und y Koordinaten eines Elements
% Input: Lame Parameter mu und lambda

% Output: Lokale Steifigkeitsmatrix

% [area,b,c]=hatGradients(x,y); % Fläche des Elements und Gradienten der Basisfkt. bestimmen
% D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0]; % Elastizität Matrix aufstellen
% % Spannungsmatrix aufstellen
% BK=[b(1) 0 b(2) 0 b(3) 0 ;
%     0 c(1) 0 c(2) 0 c(3);
%     c(1) b(1) c(2) b(2) c(3) b(3)];
% KK=BK'*D*BK*area; % Lokale Steifigkeitsmatrix bestimmen


[area,phi_jacobi]=hatGradients(x,y,d_phihat); % Fläche des Elements und Jacobi-Matix der Basisfkt. bestimmen
D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0]; % Elastizität Matrix aufstellen
% Verzerrungsmatrix aufstellen %TODO: kein hard-coding?
BK=[phi_jacobi(1,1),0,              phi_jacobi(2,1),0,              phi_jacobi(3,1),0;
   0,              phi_jacobi(1,2),0,              phi_jacobi(2,2),0,              phi_jacobi(3,2);
   phi_jacobi(1,2),phi_jacobi(1,1),phi_jacobi(2,2),phi_jacobi(2,1),phi_jacobi(3,2),phi_jacobi(3,1)];
%fuer konstante Gradienten TODO: fuer nicht konstante anpassen
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

% function [area,b,c] = hatGradients(x,y)
% % Input: x und y Koordinaten eines Elements
% 
% % Output: ...
% area=polyarea(x,y); % Fläche mittels polyarea bestimmen
% b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area; % Gradient bestimmen
% c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area; % Gradient bestimmen
% end

function [area,phi_jacobi] = hatGradients(x,y,d_phihat)
% Input: x und y Koordinaten eines Elements

% Output: ...
area=polyarea(x,y); % Fläche mittels polyarea bestimmen

%Affin lineare Abbildung
[B_affmap,d_affmap] = aff_map(x,y);
%detB_affmap = abs(det(B_affmap));
InvB_affmap = B_affmap\eye(size(B_affmap));

phi_jacobi=zeros(length(x),2);
for i=1:length(x)
    %TODO: hier erstmal funktionsauswertungen egal, da konstante gradienten
    phi_jacobi(i,:)=[d_phihat{1,i}(x,y), d_phihat{2,i}(x,y)]*InvB_affmap;
end

end

%% Affine mapping function
function [B_affmap,d_affmap] = aff_map(x,y)
a1 = [x(1);y(1)];
a2 = [x(2);y(2)];
a3 = [x(3);y(3)];

d_affmap = a1;
B_affmap = [a2-a1 , a3-a1];
end