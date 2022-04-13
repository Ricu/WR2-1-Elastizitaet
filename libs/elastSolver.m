function [U,V] = elastSolver(grid,E,nu,force,gD)
vert = grid.vert; tri = grid.tri; dirichlet = grid.dirichlet;
tri = [tri,ones(length(tri),1)];
e = 0;
ind = [2*find(dirichlet)-1 ; 2*find(dirichlet)];
dirichlet = false(2*length(vert),1);
dirichlet(ind) = true;
[mu,lambda]=enu2lame(E,nu);
[K,~,F]=elastAssemble(vert',e,tri',mu,lambda,force);


d = zeros(length(F),1); % nodal displacement vector
d(dirichlet) = gD(vert(dirichlet));
F = F(~dirichlet)- K(~dirichlet,dirichlet)*d(dirichlet); % modify load for BC
K = K(~dirichlet,~dirichlet); % modify stiffness for BC
d(~dirichlet) = K\F; % solve for free DoFs
U = d(1:2:end); 
V = d(2:2:end);
end

function [K,M,F,K_full] = elastAssemble(p,e,t,lambda,mu,force,dirichlet,gD)
% Input: p als matrix mit allen Knoten
% Input: e als matrix mit allen Kanten
% Input: t als matrix mit allen Verbindungen der Trainagulierung

% Output: K als global assemblierte Steifigkeitsmatrix
% Output: M als global assemblierte Massenmatrix
% Output: F als global assemblierter load vector
ndof=2*size(p,2); % absolute Anzahl der Freiheitsgrade entspricht ...
% der Anzahl an Knoten*2
K=sparse(ndof,ndof); % initialisiere Steifigkeitsmatrix
M=sparse(ndof,ndof); % initialisiere Massenmatrix
F=zeros(ndof,1); % initialisiere load vector
dofs=zeros(6,1); % initialisiere Anzahl Freiheitsgrade je Element
for i=1:size(t,2) % Elementweises Vorgehen
    nodes=t(1:3,i); % Elementknoten
    vert=p(1,nodes); y=p(2,nodes); % Koordinaten der Knoten
    
    f=force(vert,y); % evaluate force at nodes
    KK=elasticStiffness(vert,y,lambda,mu); % Element Steifigkeitsmatrix
    MK=elasticMass(vert,y); % Element Massenmatrix
    fK=[f(1,1) f(2,1) f(1,2) f(2,2) f(1,3) f(2,3)]'; % nodal force ...
    % values aufstellen
    FK=MK*fK; % Element load vector
    
    % Assembliere durch addieren der Element Matrixen an die richtigen
    % Stellen der globalen Matrixen
    % The two displacement components in node number i is mapped onto
    % vector entries d(2i-1) and d(2i)
    dofs(2:2:end)=2*nodes; % element degrees of freedom
    dofs(1:2:end)=2*nodes-1;
    K(dofs,dofs)=K(dofs,dofs)+KK; % addieren zur Steifigkietsmatrix
    M(dofs,dofs)=M(dofs,dofs)+MK; % addieren zur Massenmatrix
    F(dofs)=F(dofs)+FK; % addieren zum load vector
end

end


function KK = elasticStiffness(vert,y,mu,lambda)
% Input: vert und y Koordinaten eines Elements
% Input: Lame Parameter mu und lambda

% Output: Element Steifigkeitsmatrix

% triangle area and gradients (b,c) of hat functions
[area,b,c]=hatGradients(vert,y);
% elastic matrix
D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0];
% strain matrix
BK=[b(1) 0 b(2) 0 b(3) 0 ;
    0 c(1) 0 c(2) 0 c(3);
    c(1) b(1) c(2) b(2) c(3) b(3)];
% element stiffness matrix
KK=BK'*D*BK*area;
end

function MK = elasticMass(vert,y)
area=polyarea(vert,y);
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
mu=E/(2*(1+nu));
lambda=E*nu/((1+nu)*(1-2*nu));
end

function f = force(vert,y)
f=[35/13*y-35/13*y.^2+10/13*vert-10/13*vert.^2;
    -25/26*(-1+2*y).*(-1+2*vert)];
end

function [area,b,c] = hatGradients(vert,y)
% Input: vert und y Koordinaten eines Elements

% Output: ...
area=polyarea(vert,y);
b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c=[vert(3)-vert(2); vert(1)-vert(3); vert(2)-vert(1)]/2/area;
end