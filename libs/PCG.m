function [x,k,alpha,beta] = PCG(A,b,invM)

%Initialisierung
x=zeros(length(b),1);
r0=b-A*x;
z0=invM(r0);
p=z0;

%Speicherreservierung fuer alpha und beta
alpha=zeros(1000,1);
beta=zeros(1000,1);

%Iteration ueber k bis Genauigkeit erreicht
k=0;
while norm(r0)/norm(x) > 10^-8
    alphak=(r0'*z0)/(p'*A*p);
    x=x+alphak*p;
    r1=r0-alphak*A*p;
    z1=invM(r1);
    betak=(z1'*r1)/(z0'*r0);
    p=z1+betak*p;
    
    r0=r1;
    z0=z1;
    
    k=k+1;
    alpha(k)=alphak;
    beta(k)=betak;
end

alpha=alpha(1:k);
beta=beta(1:k);

end

