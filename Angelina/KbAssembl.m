function [K,b] = KbAssembl(f,daten,gd,tri,knoten,numKnlok,numKn,numEl,D,I)

%P1-Basisfunktionen
[phi_hat_p1,Dphi_hat_p1]=P1_Basisfkt;

K=sparse(numKn,numKn);
b=zeros(numKn,1);

for t=1:numEl
    T=tri(t,:);
    Tx=knoten(T,:);
    
    [B,d] = affinLinAbb(Tx');
    Binv=inv(B);
    absDetB=abs(det(B));
    
    K_T=zeros(numKnlok);
    b_T=zeros(numKnlok,1);
    
    for i=1:numKnlok
        for j=1:numKnlok
            %Elementsteifigkeitsmatrix
            for k=1:length(daten.quadratur_P2.gewichte)
                wk=daten.quadratur_P2.gewichte(k);
                xk=daten.quadratur_P2.knoten(k,:);
                K_T(i,j)=K_T(i,j)+wk*absDetB*dot(Dphi_hat_p1{i}(xk)*Binv,Dphi_hat_p1{j}(xk)*Binv);
            end
        end
        %Elementlastvektor
        for k=1:length(daten.quadratur_P2.gewichte)
            wk=daten.quadratur_P2.gewichte(k);
            xk=daten.quadratur_P2.knoten(k,:);
            xk_T=B*xk'+d;
            b_T(i)=b_T(i)+wk*absDetB*dot(phi_hat_p1{i}(xk),f(xk_T));
        end
    end
    
    %  4. Assemblierung: Globale Matrizen und Vektoren
    K(T,T)=K(T,T)+K_T;
    b(T)=b(T)+b_T;
end

%5. Randbedingungen
b(D)=gd(knoten(D,:));
b(I)=b(I)-K(I,D)*b(D);
K(D,:)=0;
K(:,D)=0;
K(D,D)=speye(nnz(D));

end




