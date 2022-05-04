function [x,tri] = genMeshSquare(N,n)
% Gitter fuer [0,1]x[0,1] erzeugen (Dreiecke).
%
% Input: Anzahl Teilgebiete: N^2
% Input: Anzahl Knoten pro Teilgebiet: (n+1)^2
%
% Output: Punkteliste x,  s x 2 (s: Anzahl Punkte)
% Output: Elementliste tri, r x 3 (r: Anzahl Dreiecke)

assert(n > 0 && N > 0)

%% Punktegitter erzeugen
linGitter = linspace(0,1,N*n+1);
[xx,yy] = meshgrid(linGitter,linGitter);
x = [xx(:),yy(:)];

%% Elemente erzeugen
% Fuer diesen Spezialfall kann man ein Gitter sehr einfach mittels
% einer Delaunay-Triangulierung erzeugen.
% Im Allgemeinen funktioniert dies jedoch nicht!
tri = delaunay(x(:,1),x(:,2));
end