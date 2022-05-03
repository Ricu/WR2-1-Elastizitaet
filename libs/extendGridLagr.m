function [vert,tri] = extendGridLagr(vert,tri,order)
% extendGridLagr erweitert das Gitter fuer P1-Elemente auf P2-Elemente
% indem die Kantenmittelpunkte hinzugefuegt werden

% Input: vert als Knotenliste
% Input: tri als Elementliste mit 3 Knoten je Element
% Input: order als gewuenschte Ordnung der Elemente

% Output: vert als erweiterte Knotenliste
% Output: tri als erweiterte Elementliste mit z.B. 6 Knoten je Element fuer
%         Order = 1

if order == 1 % Fuer order = 1 gibt es nicht zu tun.
    return
end

nElements = size(tri,1);            % Anzahl Elemente
maxIndOrig = size(vert,1);          % Hoechster urspruenglicher Knotenindex
nAdd = sum(3:order+1);              % Anzahl an Knoten die pro Element hinzugefuegt werden
vertNew = zeros(nAdd*nElements,2);  % Neue Knotenliste (inkl. moeglicher Duplikate)
counter = 1;                        % Zaehler

% Idee: berechne fuer jedes Element die Knotenmittelpunkte und eliminiere
% anschliessend die doppelt berechnenten Knoten. Mittels eines dabei
% entstehenden Abbildungsvektors lassen sich die Elementknoten einfach
% zuordnen. Berechne dies alles fuer die neuen Knoten und fuege die neue
% Information an die alten an.
for i = 1:nElements
    x = vert(tri(i,:),1);               % x-Koordinaten der Elementknoten
    y = vert(tri(i,:),2);               % y-Koordinaten der Elementknoten
    window = counter:counter+nAdd-1;    % aktuell betrachtete Indizes
    % Berechne die Knotenmittelpunkte in x und y Komponente mittels circshift
    vertNew(window,1) = 1/2*(x+circshift(x,1)); % x-Komponente
    vertNew(window,2) = 1/2*(y+circshift(y,1)); % y-Komponente
    counter = counter + nAdd;                   % Aktualisiere den Zaehler
end

% Entferne die Knoten, welche mehrfach erstellt wurden. 
% ic bildet die Knotenindizes vor der Eliminierung auf jene danach ab.
% Dies kann dann im Folgenden zur Markierung der Knoten auf die enzelnen
% Elemente genutzt werden. Praktischerweise wurden die urspruenglichen
% Knoten in der richtigen Reihenfolge erstellt 
[vertNew,~,ic] = unique(vertNew,"rows","stable"); 

% Bringe den Abbildungsvektor in die Form fuer die Elementliste
triNew = reshape(ic, [nAdd,size(tri,1)])'; 
% Korrigiere den Index um die bereits bestehenden Knoten
triNew = triNew + maxIndOrig;
% Fuege die neuen Knoten an die alten an
vert = [vert; vertNew];
% Fuege die neue Elementmarkierung an die alte an
tri=[tri,triNew];
end
