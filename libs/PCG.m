function [x,k,alpha,beta] = PCG(A,b,invM)
% Input: Systemmatrix A
% Input: Rechte Seite b
% Input: Vorkonditionierer invM

% Output: Loesung des LGS x
% Output: Anzahl Iterationen k
% Output: Parameter der Iteration: alpha und beta

%% Initialisierung
x=zeros(length(b),1); % Loesungsvektor initialisieren
r0=b-A*x;   % Residuum
z0=invM*r0; % Vorkonditioniertes Residuum     
p=z0;       % Abstiegsrichtung

%Speicherreservierung fuer alpha und beta
alpha=zeros(1000,1);
beta=zeros(1000,1);

%% Iteration bis geforderte Genauigkeit erreicht
k=0; % Anzahl Iterationen
while norm(r0)/norm(x) > 10^-8 % Abbruchbedingung
    alphak=(r0'*z0)/(p'*A*p);
    x=x+alphak*p;       % Loesungsvektor
    r1=r0-alphak*A*p;   % neues Residuum
    z1=invM*r1;         % vorkonditioniertes neues Residuum
    betak=(z1'*r1)/(z0'*r0);
    p=z1+betak*p;       % Abstiegsrichtung
    
    % Aktualisierung der Werte fuer naechste Iteration
    r0=r1;
    z0=z1;
    
    k=k+1;
    alpha(k)=alphak;
    beta(k)=betak;
end

% Entferne nicht benoetigten Speicherplatz fuer alpha und beta
alpha=alpha(1:k);
beta=beta(1:k);

end

