function [vert,tri] = extendGridLagr(vert,tri,order)
%EXTENDGRIDLAGR Summary of this function goes here
%   Detailed explanation goes here
% Annahme, dass hier gerade order 2 verwendet wird
nElements = size(tri,1);
vertCounter = size(vert,1);
nAdd = 3; % Hardcode fuer order = 2; Anzahl an Knoten die Pro element hinzugefügt werden
vertNew = zeros(nAdd*nElements,2); % Neue Knotenliste (inkl. moeglicher Duplikate)
counter = 1;

for i = 1:nElements
    x = vert(tri(i,:),1);
    y = vert(tri(i,:),2);
    window = counter:counter+nAdd-1;
    vertNew(window,1) = 1/2*(x+circshift(x,1)); %circshift verschiebt element um 1
    vertNew(window,2) = 1/2*(y+circshift(y,1)); %circshift verschiebt element um 1
    counter = counter + nAdd;
end

[vertNew,~,ic] = unique(vertNew,"rows","stable"); %Entferne die Knoten, welche mehrfach erstellet wurden
% ic ist nun eine Abbildung von mehrfachen Knoten auf die einzigartigen
% dies kann dann im Folgenden zur Markierung der Knoten auf die einzelnen
% Elemente genutzt werden. Praktischerweise wurden ja die urspruenglichen
% Knoten in der richtigen Reihenfolge erstellt (in dem loop)
triNew1 = reshape(ic, [nAdd,size(tri,1)])'; %markiere Knoten
triNew = triNew1 + vertCounter; %Korrigiere um die bereits bestehenden Knoten
vert = [vert; vertNew];
tri = [tri, triNew];
end

