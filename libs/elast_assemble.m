function [K,M,b] = elast_assemble(tri,vert,order,f)
%Assemblierungsroutine
% Input: 
% tri: 
% vert:
% order:
% f:

% Output:
% K: Steifigkeitsmatrix fuer das gegebene Element
% M: Massenmatrix fuer das gegebene Element
% b: Lastvektor fuer das gegebene Element

% Number of x and base functions
numElements = size(tri,1);
numBaseFun = (order+2)*(order+1)/2;
nBF2 = numBaseFun^2; % Anzahl Elemente der lokalen Steifigkeitsmatrix
gmres()
% Prepare matrices for sparse()
K_val = zeros(nBF2*numElements,1); 
M_val = K_val;
iIndex = K_val; 
jIndex = K_val;
b = zeros(size(vert,1),1);

% Load function and quadrature data
[phi,d_phi] = baseFun(order);
if order == 1
    quad_low = load("quad_formeln.mat").quadratur_P1;
    quad_high = load("quad_formeln.mat").quadratur_P5;
elseif order == 2
    quad_low = load("quad_formeln.mat").quadratur_P2;
    quad_high = load("quad_formeln.mat").quadratur_P5;
end

for i = 1:size(tri,1)
    [B,d]=aff_map(vert,tri(i,:));
    [K_T,M_T,b_T] = getMatrices(B,d,f,phi,d_phi,quad_low,quad_high);
    
    % Assignments for sparse()
    iIndex((i-1)*nBF2+1:i*nBF2) = reshape(repmat(tri(i,:),numBaseFun,1),nBF2,1);
    jIndex((i-1)*nBF2+1:i*nBF2) = repmat(tri(i,:),1,numBaseFun); 
    K_val((i-1)*nBF2+1:i*nBF2) = reshape(K_T,nBF2,1);
    M_val((i-1)*nBF2+1:i*nBF2) = reshape(M_T,nBF2,1);
    b(tri(i,:)) = b(tri(i,:)) + b_T;
end
n = size(vert,1);
K = sparse(iIndex,jIndex,K_val,n,n);
M = sparse(iIndex,jIndex,M_val,n,n);

%TODO Randbedingung hier
% Wenn in der Funktion: brauchen argumente:
% dirichlet markiert Dirichletrand
% inner markiert Innere Knoten
% gDirichlet ist Dirichletrand funktion
% 
% b(dirichlet)=gDirichlet(knoten(dirichlet,:));
% b(inner)=b(inner)-K(inner,dirichlet)*b(dirichlet);
% K(dirichlet,:)=0;
% K(:,dirichlet)=0;
% K(dirichlet,dirichlet)=speye(nnz(dirichlet));

end

%% Get matrices
function [K,M,b] = getMatrices(B,d,f,phi,d_phi,quad_low,quad_high)

numBaseFunc = length(phi);
K = zeros(numBaseFunc);
M = zeros(numBaseFunc);
b = zeros(numBaseFunc,1);

detb = abs(det(B));
invb = B\eye(size(B));

for i = 1:numBaseFunc
    for j = 1:numBaseFunc
        %% Matrix K_T
        x = quad_low.knoten(:,1);
        y = quad_low.knoten(:,2);
        %TODO Effizienz vs. Intuition
        temp = dot([d_phi{1,i}(x,y), d_phi{2,i}(x,y)]*invb,...
                   [d_phi{1,j}(x,y), d_phi{2,j}(x,y)]*invb,2);
        K(i,j) = detb* sum(quad_low.gewichte .* temp);
    end  
    
    %% Vektor b_T
    v_ref = quad_high.knoten';
    v = B*v_ref +d;
    %TODO Effizienz vs. Intuition
    temp = f(v(1,:),v(2,:)).*phi{i}(v_ref(1,:),v_ref(2,:));
    b(i) = detb* sum(quad_high.gewichte'.*temp);
end
end

%% Affine mapping function
function [B,d] = aff_map(coords,x)
a1 = coords(x(1),1:2)';
a2 = coords(x(2),1:2)';
a3 = coords(x(3),1:2)';

d = a1;
B = [a2-a1 , a3-a1];
end

function ElasticSolver()
    g=Rectg(0,0,1,1);
    [p,e,t]=initmesh(g,"hmax",0.1);
    E=1; nu=0.3;
    [mu,lambda]=Enu2Lame(E,nu);
    [K,M,F]=ElasticAssembler(p,e,t,mu,lambda,@Force);
    bdry=unique([e(1,:) e(2,:)]);
    fixed=[2*bdry-1 2*bdry];
    values=zeros(length(fixed),1);
    ndof=length(F);
    free=setdiff([1:ndof],fixed);
    F=F(free)-K(free,fixed)*values;
    K=K(free,free);
    d=zeros(ndof,1);
    d(free)=K\F;
    d(fixed)=values;
    U=d(1:2:end); V=d(2:2:end);
    figure(1), pdesurf(p,t,U), title("(u_h)_1")
    figure(2), pdesurf(p,t,V), title("(u_h)_2")
end

function [K,M,F] = ElasticAssembler(p,e,t,lambda,mu,force)
    ndof=2*size(p,2); % total number of degrees of freedom
    K=sparse(ndof,ndof); % allocate stiffness matrix
    M=sparse(ndof,ndof); % allocate mass matrix
    F=zeros(ndof,1); % allocate load vector
    dofs=zeros(6,1); % allocate element degrees of freedom
    
    for i=1:size(t,2) % assemly loop over elements
    nodes=t(1:3,i); % element nodes
    x=p(1,nodes); y=p(2,nodes); % node coordinates
    dofs(2:2:end)=2*nodes; % element degrees of freedom
    dofs(1:2:end)=2*nodes-1;
    f=force(x,y); % evaluate force at nodes
    KK=ElasticStiffness(x,y,lambda,mu); % element stiffness
    MK=ElasticMass(x,y); % element mass
    fK=[f(1,1) f(2,1) f(1,2) f(2,2) f(1,3) f(2,3)]’;
    FK=MK*fK; % element load
    K(dofs,dofs)=K(dofs,dofs)+KK; % add to stiffness matrix
    M(dofs,dofs)=M(dofs,dofs)+MK; % add to mass matrix
    F(dofs)=F(dofs)+FK; % add to load vector
    end
    
    bdry=unique([e(1,:) e(2,:)]); % boundary nodes
    fixed=[2*bdry-1 2*bdry]; % boundary degrees of freedom, DoFs
    values=zeros(length(fixed),1); % zero boundary values
    ndof=length(F); % total number of DoFs
    free=setdiff([1:ndof],fixed); % free DoFs
    F=F(free)-K(free,fixed)*values; % modify load for BC
    K=K(free,free); % modify stiffness for BC
    d=zeros(ndof,1); % nodal displacement vector
    d(free)=K\F; % solve for free DoFs
    d(fixed)=values; % insert known DoFs
end


function KK = ElasticStiffness(x,y,mu,lambda)
%     % triangle area and gradients (b,c) of hat functions
%     [area,b,c]=HatGradients(x,y);
%     % elastic matrix
%     D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0];
%     % strain matrix
%     BK=[b(1) 0 b(2) 0 b(3) 0 ;
%     0 c(1) 0 c(2) 0 c(3);
%     c(1) b(1) c(2) b(2) c(3) b(3)];
%     % element stiffness matrix
%     KK=BK'*D*BK*area;

for t=1:numEl %%TODO Variable numEl
    [B_aL,d_aL] = affinLinAbb(x,y);
    
    B_aLinv=inv(B_aL);
    absDetB=abs(det(B_aL));
    
%     K_T=zeros(numKnlok);
%     b_T=zeros(numKnlok,1);
%     
%     for i=1:numKnlok
%         for j=1:numKnlok
%             %Elementsteifigkeitsmatrix
%             for k=1:length(daten.quadratur_P2.gewichte)
%                 wk=daten.quadratur_P2.gewichte(k);
%                 xk=daten.quadratur_P2.knoten(k,:);
%                 K_T(i,j)=K_T(i,j)+wk*absDetB*dot(Dphi_hat_p1{i}(xk)*Binv,Dphi_hat_p1{j}(xk)*Binv);
%             end
%         end
%         %Elementlastvektor
%         for k=1:length(daten.quadratur_P2.gewichte)
%             wk=daten.quadratur_P2.gewichte(k);
%             xk=daten.quadratur_P2.knoten(k,:);
%             xk_T=B*xk'+d;
%             b_T(i)=b_T(i)+wk*absDetB*dot(phi_hat_p1{i}(xk),f(xk_T));
%         end
%     end
    
    
end



D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0];
KK=BK'*D*BK*area;

end

function MK = ElasticMass(x,y)
    area=polyarea(x,y);
    MK=[2 0 1 0 1 0;
    0 2 0 1 0 1;
    1 0 2 0 1 0;
    0 1 0 2 0 1;
    1 0 1 0 2 0;
    0 1 0 1 0 2]*area/12;
end

function [mu,lambda] = Enu2Lame(E,nu)
    mu=E/(2*(1+nu));
    lambda=E*nu/((1+nu)*(1-2*nu));
end

function f = Force(x,y)
    f=[35/13*y-35/13*y.^2+10/13*x-10/13*x.^2;
    -25/26*(-1+2*y).*(-1+2*x)];
end

function [B_Al,d_Al] = affinLinAbb(x,y)

B_Al=zeros(2);
B_Al(1,1)=x(2)-x(1);
B_Al(2,1)=y(2)-y(1);
B_Al(1,2)=x(3)-x(1);
B_Al(2,2)=y(3)-y(1);

d_Al=zeros(2,1);
d_Al(1)=x(1);
d_Al(2)=y(1);

end
