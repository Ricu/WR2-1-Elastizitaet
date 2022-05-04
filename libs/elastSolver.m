function [U,V,K,F] = elastSolver(grid,E,nu,f,gD,order)
% Input: grid Structure mit allen Gitterkomponenten
% Input: E der Young modulus
% Input: nu die Poisson ratio
% Input: f Volumenkraft
% Input: gD Funktion des Dirichletrandes
% Input: order Ordnung der Elemente

% Output: U Deformation in x1-Richtung
% Output: V Deformation in x2-Richtung
% Output: K Steifigkeitsmatrix mit eliminertem Dirichletrand
% Output: F Lastvektor mit eliminiertem Dirichletrand

vert = grid.vert; tri = grid.tri; dirichlet = grid.dirichlet; % Gitterdaten extrahieren
dirichlet2 = repelem(dirichlet,2);          % Erweitere den Dirichletrand um die 2. Komponente

[mu,lambda]=enu2lame(E,nu);                 % Materialparameter extrahieren
[K,~,F] = elastAssemble(vert',tri',mu,lambda,f,dirichlet2,gD,order); % Assemblierungsroutine aufrufen

d = zeros(length(F),1);                     % Loesungsvektor initialisieren
d(dirichlet2) = gD(vert(dirichlet2));       % Dirichletwerte eintrage
d(~dirichlet2) = K\F;                       % Galerkin System loesen
U = d(1:2:end);                             % Loesung in x_1 Richtung extrahieren
V = d(2:2:end);                             % Loesung in x_2 Richtung extrahieren
end

function [K,M,F] = elastAssemble(p,t,lambda,mu,f,dirichlet,gD,order)
% Input: p als matrix mit allen Knoten
% Input: e als matrix mit allen Kanten
% Input: t als matrix mit allen Verbindungen der Triangulierung

% Output: K als global assemblierte Steifigkeitsmatrix
% Output: M als global assemblierte Massenmatrix
% Output: F als global assemblierter Lastvektor

%% Lade Basisfunktionen und Quadraturformeln
[phihat,d_phihat] = baseFun(order);
if order == 1
    quad_low = load("quad_formeln.mat").quadratur_P1;
    quad_high = load("quad_formeln.mat").quadratur_P5;
elseif order == 2
    quad_low = load("quad_formeln.mat").quadratur_P2;
    quad_high = load("quad_formeln.mat").quadratur_P5;
end

%% Matrizen vorbereiten, in die die Index-Wertepaare gespeichert werden
nBaseFun = 2*length(phihat);        % Anzahl Basisfunktionen
nElem = length(t);                  % Anzahl Elemente global
nValLoc = nBaseFun^2;               % Anzahl Eintraege der lokalen Steifigkeitsmatrix
K_val = zeros(nValLoc*nElem,1);     % Liste fuer K-Werte
M_val = K_val;                      % Liste fuer M-Werte
iIndex = K_val;                     % Liste fuer i-Indizes
jIndex = K_val;                     % Liste fuer j-Indizes
n_nodes = 2*size(p,2);              % Anzahl Freiheitsgrade entspricht Anzahl Knoten*2
F=zeros(n_nodes,1);                 % Initialisiere Lastvektor

%% Assemblierungsroutine
for i=1:size(t,2)                               % Ueber die Elemente iterieren
    node_ind=t(1:nBaseFun/2,i);                 % Die zum aktuellen Element gehoerenden physikalischen Knoten Indizes
    x=p(1,node_ind); y=p(2,node_ind);           % Koordinaten der Knoten
       
    [B_affmap,d_affmap] = aff_map(x,y);         % Affine Abbildung aufstellen
    detB_affmap = abs(det(B_affmap));           % Determinante der affinen Abb.
    InvB_affmap = B_affmap\eye(size(B_affmap)); % Inverse der Abbildungsmatrix
    
    % Lokale Matrizen aufstellen
    KK = elasticStiffness(lambda,mu,d_phihat,quad_low,nBaseFun,InvB_affmap,detB_affmap);
    [MK,FK] = elasticMassLoad(phihat,f,B_affmap,d_affmap,detB_affmap,nBaseFun,quad_high);
    
    % Indizes plus zugehoerige Werte an die richtigen Stellen schreiben
    node_ind2 = 2*repelem(node_ind,2);           % Um die 2. Komponente erweitern
    node_ind2(1:2:end) = node_ind2(1:2:end) - 1; % Korrigiere Indizes
   
    iIndex((i-1)*nValLoc+1:i*nValLoc) = reshape(repmat(node_ind2',nBaseFun,1),nValLoc,1); % i-Indizes speichern
    jIndex((i-1)*nValLoc+1:i*nValLoc) = repmat(node_ind2,nBaseFun,1); % j-Indizes speichern
    K_val((i-1)*nValLoc+1:i*nValLoc) = reshape(KK,nValLoc,1);   % Werte für die Steifigkeitsmatrix speichern
    M_val((i-1)*nValLoc+1:i*nValLoc) = reshape(MK,nValLoc,1);   % Werte für die Massenmatrix speichern
    F(node_ind2) = F(node_ind2) + FK;                           % Der Lastvektor kann ohne Umwege aktualisiert werden
end
% Globale Matrizen aufstellen
K = sparse(iIndex,jIndex,K_val,n_nodes,n_nodes); % Anhand der zuvor erstellten Listen die Steifigkeitsmatrix erstellen
M = sparse(iIndex,jIndex,M_val,n_nodes,n_nodes); % Anhand der zuvor erstellten Listen die Massenmatrix erstellen


double_points = repelem(p,2,1);                  % Knotenliste um 2. Komponente erweitern
% Dirichletknoten eliminieren
F = F(~dirichlet)- K(~dirichlet,dirichlet)*gD(double_points(dirichlet)); 
K = K(~dirichlet,~dirichlet);
end


function KK = elasticStiffness(mu,lambda,d_phihat,quad_low,nBaseFun,InvB_affmap,detB_affmap)
% Input: Lame Parameter mu und lambda
% Input: d_phihat partiellen Ableitungen der Basisfunktionen auf dem
%        Referenzelement
% Input: quad_low Quadraturformeln für ... ?
% Input: nBaseFun Basisfunktionen entsprechender Ordnung
% Input: InvB_affmap inverse Matrix der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement 
% Input: detB_affmap Determinante der Matrix der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement

% Output: lokale Steifigkeitsmatrix

% Aufstellen der Elastizitaetsmatrix
D = mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0];

% Aufstellen der Jacobi-Matrix der Basisfunktionen des physikalischen Elements mit Transformation
J = cell(nBaseFun,2);
for i = 1:nBaseFun
    for j = 1:2
        J{i,j} = @(x,y) [d_phihat{1,i}(x,y),d_phihat{2,i}(x,y)]*InvB_affmap(:,j);
    end
end

% Aufstellen der Elementverzerrungsmatrix
if nBaseFun == 6 % P1 Elemente
    BKhat=@(x,y)[J{1,1}(x,y),0,          J{2,1}(x,y),0,          J{3,1}(x,y),0;
                 0,          J{1,2}(x,y),0,          J{2,2}(x,y),0,          J{3,2}(x,y);
                 J{1,2}(x,y),J{1,1}(x,y),J{2,2}(x,y),J{2,1}(x,y),J{3,2}(x,y),J{3,1}(x,y)];
elseif nBaseFun == 12 % P2 Elemente
    BKhat=@(x,y)[J{1,1}(x,y),0,          J{2,1}(x,y),0,          J{3,1}(x,y),0,          J{4,1}(x,y),0,          J{5,1}(x,y),0,          J{6,1}(x,y),0;
                 0,          J{1,2}(x,y),0,          J{2,2}(x,y),0,          J{3,2}(x,y),0,          J{4,2}(x,y),0,          J{5,2}(x,y),0,          J{6,2}(x,y);
                 J{1,2}(x,y),J{1,1}(x,y),J{2,2}(x,y),J{2,1}(x,y),J{3,2}(x,y),J{3,1}(x,y),J{4,2}(x,y),J{4,1}(x,y),J{5,2}(x,y),J{5,1}(x,y),J{6,2}(x,y),J{6,1}(x,y)];
end

%% Aufstellen der Elementstifigkeismatrix: Loesen des Integrals mit Quadraturformel
% Knoten der Quadraturformel
xhat_quad = quad_low.knoten(:,1);  yhat_quad = quad_low.knoten(:,2);

KK=zeros(nBaseFun);
for i=1:length(xhat_quad)
    % Ausdruck an Quadraturknoten auswerten
    temp=BKhat(xhat_quad(i),yhat_quad(i))'*D*BKhat(xhat_quad(i),yhat_quad(i));
    KK = KK + detB_affmap* quad_low.gewichte(i) .* temp;
end
end

function [MK,FK]=elasticMassLoad(phihat,f,B_affmap,d_affmap,detB_affmap,nBaseFun,quad_high)
% Input: phihat Basisfunktionen auf dem Referenzelement
% Input: f Volumenkraft
% Input: B_affmap Matrix der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement
% Input: d_affmap Vektor der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement
% Input: detB_affmap Determinante der Matrix der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement
% Input: nBaseFun Basisfunktionen der entsprechenden Ordnung
% Input: quad_high Quadraturformeln fuer... ?

% Output: lokaler Lastvektor
% Output: lokale Massenmatrix

% Aufstellen der Matrix, die die Basisfunktionen des Referenzelements enthaelt
if nBaseFun == 6 % P1 Elemente
    phihat_mat=@(x,y)[phihat{1}(x,y),0;
                      0,phihat{1}(x,y);
                      phihat{2}(x,y),0;
                      0,phihat{2}(x,y);
                      phihat{3}(x,y),0;
                      0,phihat{3}(x,y)];
elseif nBaseFun == 12 % P2 Elemente
    phihat_mat=@(x,y)[phihat{1}(x,y),0;
                      0,phihat{1}(x,y);
                      phihat{2}(x,y),0;
                      0,phihat{2}(x,y);
                      phihat{3}(x,y),0;
                      0,phihat{3}(x,y);
                      phihat{4}(x,y),0;
                      0,phihat{4}(x,y);
                      phihat{5}(x,y),0;
                      0,phihat{5}(x,y);
                      phihat{6}(x,y),0;
                      0,phihat{6}(x,y);];
end

% Knoten der Quadraturformel
xhat_quad = quad_high.knoten(:,1);  % x-Komponente Quadraturknoten
yhat_quad = quad_high.knoten(:,2);  % y-Komponente Quadraturknoten
nVertQuad=length(xhat_quad);        % Anzahl der Quadraturknoten

% Aufstellen der Elementmassenmatrix: Loesen des Integrals mit geeigneter Quadraturformel
MK=zeros(nBaseFun);
for i=1:length(xhat_quad)
    % Ausdruck an Quadraturknoten auswerten
    temp=phihat_mat(xhat_quad(i),yhat_quad(i))*phihat_mat(xhat_quad(i),yhat_quad(i))';
    MK = MK + detB_affmap* quad_high.gewichte(i) .* temp;
end

% Transformation der Quadraturknoten auf das physikalische Element
x_quad=zeros(nVertQuad,1); y_quad=zeros(nVertQuad,1);
for i=1:nVertQuad
    x_quad(i)=B_affmap(1,:)*[xhat_quad(i);yhat_quad(i)]+d_affmap(1);
    y_quad(i)=B_affmap(2,:)*[xhat_quad(i);yhat_quad(i)]+d_affmap(2);
end

% Aufstellen des Elementlastvektors: Loesen des Integrals mit geeigneter Quadraturformel
FK=zeros(nBaseFun,1);
for i=1:length(xhat_quad)
    % Ausdruck an Quadraturknoten auswerten
    temp=phihat_mat(xhat_quad(i),yhat_quad(i))*f(x_quad(i),y_quad(i));
    FK = FK + detB_affmap* quad_high.gewichte(i) .* temp;
end
end

function [mu,lambda] = enu2lame(E,nu)
% Input: E ist der Young modulus
% Input: nu ist die Poisson ratio

% Output: Lame Parameter mu und lambda
mu=E/(2*(1+nu));                    % mu berechnen
lambda=E*nu/((1+nu)*(1-2*nu));      % Lambda berechnen
end

%% Affin lineare Funktion
function [B_affmap,d_affmap] = aff_map(x,y)
%Input: x- und y-Koordinaten der Knoten eines Elements der Triangulierung

%Output: B_affmap Matrix der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement
%Output: d_affmap Vektor der affin linearen Abbildung vom physikalischen 
%        Element auf das Referenzelement

a1 = [x(1);y(1)]; a2 = [x(2);y(2)]; a3 = [x(3);y(3)]; % Knoten aufstellen
d_affmap = a1; B_affmap = [a2-a1 , a3-a1]; % Abbildung berechnen
end