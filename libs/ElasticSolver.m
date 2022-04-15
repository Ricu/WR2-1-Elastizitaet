function ElasticSolver()
%     g=Rectg(0,0,1,1);
%     [p,e,t]=initmesh(g,"hmax",0.1);
    % Gittergenerierung und Erstellung von Knoten- und Elementlisten
    [x,tri] = genMeshSquare(1,10);
    tri = [tri,ones(length(tri),1)];
    e = 0; % Kanten Matrix
    % Paramater aus der Aufgabenstellung initialisieren
    E=1; nu=0.3;
    xLim = [0,1];
    yLim = [0,1];
    % Definiere den Dirichlet-Rand
    dirichlet = or(ismember(x(:,1),xLim), ismember(x(:,2),yLim));
    ind = [2*find(dirichlet)-1 ; 2*find(dirichlet)];
    dirichlet = false(2*length(x),1);
    dirichlet(ind) = true;
    % Berechne die Lame Parameter
    [mu,lambda]=Enu2Lame(E,nu);
    % Brechne die globalen Steifigkeits-,Massenmatrix und den load vector
    [K,M,F]=ElasticAssembler(x',e,tri',mu,lambda,@Force);

    % Integriere den homogenen Dirichletrand (Haben wir schon)
    
%     bdry=unique([e(1,:) e(2,:)]); % boundary nodes
%     fixed=[2*bdry-1 2*bdry]; % boundary degrees of freedom, DoFs
%     values=zeros(length(fixed),1); % zero boundary values
%     ndof=length(F); % total number of DoFs
%     free=setdiff([1:ndof],fixed); % free DoFs
    d=zeros(length(F),1); % nodal displacement vector
    F=F(~dirichlet)-K(~dirichlet,dirichlet)*d(dirichlet); % modify load for BC
    K=K(~dirichlet,~dirichlet); % modify stiffness for BC
%     d=zeros(ndof,1); % nodal displacement vector
    d(~dirichlet)=K\F; % solve for free DoFs
%     d(dirichlet)=values; % insert known DoFs
    U=d(1:2:end); V=d(2:2:end);
    figure()
    subplot(1,2,1), trisurf(tri,x(:,1),x(:,2),U), title("(u_h)_1")
    subplot(1,2,2), trisurf(tri,x(:,1),x(:,2),V), title("(u_h)_2")
end

function [K,M,F] = ElasticAssembler(p,e,t,lambda,mu,force)
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
        x=p(1,nodes); y=p(2,nodes); % Koordinaten der Knoten

        f=force(x,y); % evaluate force at nodes
        KK=ElasticStiffness(x,y,lambda,mu); % Element Steifigkeitsmatrix
        MK=ElasticMass(x,y); % Element Massenmatrix
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


% Um sicher zu stellen, dass die Steifigkeitsmatrix invertierbar ist,
% benötigen wir Randbedingungen.

%Beispielcode für enforced homogenen Dirichlet-Rand:
    
%     bdry=unique([e(1,:) e(2,:)]); % boundary nodes
%     fixed=[2*bdry-1 2*bdry]; % boundary degrees of freedom, DoFs
%     values=zeros(length(fixed),1); % zero boundary values
%     ndof=length(F); % total number of DoFs
%     free=setdiff([1:ndof],fixed); % free DoFs
%     F=F(free)-K(free,fixed)*values; % modify load for BC
%     K=K(free,free); % modify stiffness for BC
%     d=zeros(ndof,1); % nodal displacement vector
%     d(free)=K\F; % solve for free DoFs
%     d(fixed)=values; % insert known DoFs
end


function KK = ElasticStiffness(x,y,mu,lambda)
% Input: x und y Koordinaten eines Elements
% Input: Lame Parameter mu und lambda

% Output: Element Steifigkeitsmatrix

    % triangle area and gradients (b,c) of hat functions
    [area,b,c]=HatGradients(x,y);
    % elastic matrix
    D=mu*[2 0 0; 0 2 0; 0 0 1]+lambda*[1 1 0; 1 1 0; 0 0 0];
    % strain matrix
    BK=[b(1) 0 b(2) 0 b(3) 0 ;
    0 c(1) 0 c(2) 0 c(3);
    c(1) b(1) c(2) b(2) c(3) b(3)];
    % element stiffness matrix
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
% Input: E ist der Young modulus
% Input: nu ist die Poisson ratio

% Output: Lame Parameter mu und lambda
    mu=E/(2*(1+nu));
    lambda=E*nu/((1+nu)*(1-2*nu));
end

function f = Force(x,y)
    f=[35/13*y-35/13*y.^2+10/13*x-10/13*x.^2;
    -25/26*(-1+2*y).*(-1+2*x)];
end

function [area,b,c] = HatGradients(x,y)
% Input: x und y Koordinaten eines Elements

% Output: ...
    area=polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end

function r = Rectg(xmin,ymin,xmax,ymax)
r=[2 xmin xmax ymin ymin 1 0;
2 xmax xmax ymin ymax 1 0;
2 xmax xmin ymax ymax 1 0;
2 xmin xmin ymax ymin 1 0]';
end