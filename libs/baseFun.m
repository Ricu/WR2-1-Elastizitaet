function [phi,d_phi] = baseFun(order)
vertices = [0,0; 1,0; 0,1]; % Eckknoten
if order == 1
    %% P1
    phi = {@(x,y) 1-x-y, @(x,y) x, @(x,y) y};   % Basisfkt.
    d_phi = {@(x,y) -1,@(x,y) 1, @(x,y) 0;      % Ableitung nach x
             @(x,y) -1,@(x,y) 0, @(x,y) 1};     % Ableitung nach y
elseif order == 2
    %% P2
    vertices = extendGridLagr(vertices,[1,2,3],2);  % Erweiterung um Knotenmittelpunkte
    numBaseFun = length(vertices);                  % Anzahl Basisfunktion
    phi = cell(1,numBaseFun);                       % Basisfunktionen cell
    d_phi = cell(2,numBaseFun);                     % Zeile 1: dx   Zeile 2:  dx
    A = zeros(numBaseFun);                          
    A(:,1) = ones(numBaseFun,1);                    % Konstante
    A(:,2) = vertices(:,1);                         % x
    A(:,3) = vertices(:,2);                         % y 
    A(:,4) = vertices(:,1).^2;                      % x^2
    A(:,5) = vertices(:,1).*vertices(:,2);          % xy
    A(:,6) = vertices(:,2).^2;                      % y^2
    for i = 1:numBaseFun
        b = zeros(numBaseFun,1);    b(i) = 1;       % Rechte Seite
        c = A\b;                                    % Koeffizienten berechnen
        phi{i} = @(x,y) c(1) + c(2)*x + c(3)*y + c(4)*x.^2 +c(5)*x.*y + c(6)*y.^2;
        d_phi{1,i}  = @(x,y) +c(2) + 2*c(4)*x + c(5)*y; % Ableitung nach x
        d_phi{2,i}  = @(x,y) +c(3) + c(5)*x + 2*c(6)*y; % Ableitung nach y
    end
end
end