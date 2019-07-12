function [T, Tk] = MTH_HelicSucess(s,s0,q,Tf,Tipo,njun)

r=size(q,1);
T=zeros(4,4,njun);

k=1;
%CÁLCULO DAS N*R MATRIZES DE TRANSFORMAÇÃO [T]
while k<=njun*r
    
    for i=1:r
        for j=1:njun
            T(:,:,k) = MTH_Helic(s(j,:),s0(j,:),q(i,j),Tipo(j));
            k=k+1;
        end
    end
end

%AGRUPAMENTO (PRODUTÓRIO) DAS MATRIZES DE TRANSFORMAÇÃO EM [Tk]
Tk(:,:,1)=T(:,:,1);
for k=1:r

Tk(:,:,k)=T(:,:,(k-1)*njun+1);

    for i=(k-1)*njun+1:(k*njun)-1
        Tk(:,:,k)=Tk(:,:,k)*T(:,:,i+1);
    end
        Tk(:,:,k)=Tk(:,:,k)*Tf;
end

% %EXTRAÇÃO DA POSIÇÃO E ORIENTAÇÃO DAS JUNTAS
% Junta=1;
% for Junta=1:njun
%     
% i=Junta;
%     for k=1:r
%         %POSIÇÃO
%         PJunta(Junta*)(k,:)=[T(1,4,i) T(2,4,i) T(3,4,i)];
%         %ORIENTAÇÃO
%         %...
%         
%         i=i+njun;
%     end
% 
% end
% 
% %EXTRAÇÃO DA POSIÇÃO E ORIENTAÇÃO DA FERRAMENTA
% %POSIÇÃO
% for k=1:size(Tk,3)
%     Px(k)=Tk(1,4,k);
%     Py(k)=Tk(2,4,k);
%     Pz(k)=Tk(3,4,k);
% end
% 
% %ORIENTAÇÃO
% 
% 
% end