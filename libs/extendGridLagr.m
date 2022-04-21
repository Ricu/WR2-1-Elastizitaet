function [vert,tri] = extendGridLagr(vert,tri,order)
%EXTENDGRIDLAGR Summary of this function goes here
%   Detailed explanation goes here

if order == 1
    return
end
nElements = size(tri,1);
vertCounter = size(vert,1);
nAdd = sum(3:order+1); %Anzahl an Knoten die Pro element hinzugefügt werden
vertNew = zeros(nAdd*nElements,2); % Neue Knotenliste (inkl. mglicher Duplikate)
counter = 1;

for i = 1:nElements
    x = vert(tri(i,:),1);
    y = vert(tri(i,:),2);
    window = counter:counter+nAdd-1;
    newValues = zeros(nAdd,2);
%     for j = 1:order-1
%         
%     end
%     scatter(newValues(:,1),newValues(:,2))
    vertNew(window,1) = 1/2*(x+circshift(x,1)); %circshift verschiebt element um 1
    vertNew(window,2) = 1/2*(y+circshift(y,1)); %circshift verschiebt element um 1
    counter = counter + nAdd;
end

[vertNew,~,ic] = unique(vertNew,"rows","stable"); %Entferne die Knoten, welche mehrfach erstellet wurden
% ic ist nun eine Abbildung von mehrfachen knoten auf die einzigartigen
% dies kann dann im folgende zur markierung der Knoten auf die enzelnen
% elemente genutzt werden. praktischerweise wurden ja die ursorünglichen
% knoten in der richtigen reihenfolge erstellt (in dem loop)
triNew = reshape(ic, [nAdd,size(tri,1)])'; %markiere knoten
triNew = triNew + vertCounter; %Korrigiere um die bereits bestehenden knoten
vert = [vert; vertNew];
tri=[tri,triNew];
end
