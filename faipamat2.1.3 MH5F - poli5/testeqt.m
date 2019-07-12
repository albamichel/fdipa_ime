function [dt,di,qxyz] = testeqt(x)

q_t(:,1)=x(1:10);
q_t(:,2)=x(11:20);
q_t(:,3)=x(21:30);

thetai=[0 0 0];

r=size(q_t,1);
d=0;
[s, s0, sf, s0f, njun] = Modelagem3J;
theta_f = 45*(pi/180);
Tf = MatTransHomog_f(sf,s0f,theta_f,d);

%CÁLCULO DAS N*R MATRIZES DE TRANSFORMAÇÃO
T=zeros(4,4,njun);
k=1;
while k<njun*r
    
    for i=1:r
        for j=1:njun
            T(:,:,k) = MatTransHomog(s(j,:),s0(j,:),q_t(i,j),d);
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

%O = [atan2(-T(10),T(11)); asin(T(9)); atan2(-T(5),T(1))];CORRIGIR
%Orientacao = O.*180/pi 
Posicao = [Tk(1,4,k);Tk(2,4,k);Tk(3,4,k)]

%PONTO INICIAL
Ti=eye(4,4);    
    for i=1:njun
        Tii = MatTransHomog(s(i,:),s0(i,:),thetai(i),d);
        Ti=Ti*Tii;
    end
Ti=Ti*Tf;

Qx(1)=Ti(1,4);
Qy(1)=Ti(2,4);
Qz(1)=Ti(3,4);

%EXTRAÇÃO DAS POSIÇÕES DA FERRAMENTA (X,Y,Z)
for k=2:size(Tk,3)
    Qx(k)=Tk(1,4,k);
    Qy(k)=Tk(2,4,k);
    Qz(k)=Tk(3,4,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %CÁLCULO DA DISTÂNCIA TOTAL PERCORRIDA PELA FERRAMENTA (dt)
dt=0;
    for k=2:size(Tk,3)
        di(k)=((Qx(k)-Qx(k-1))^2+(Qy(k)-Qy(k-1))^2+(Qz(k)-Qz(k-1))^2)^0.5;
    end
dt=sum(di);
di=di(:);

qxyz=[Qx' Qy' Qz'];

%SAÍDA GRÁFICA

plot(Qx,Qy,'.-')
hold on
% quiver3(X,Y,Z,ux,uy,uz,0.3)
% quiver3(X,Y,Z,vx,vy,vz,0.3)
% quiver3(X,Y,Z,wx,wy,wz,0.3)
axis([0 15 -10 10])
xlabel('X'); ylabel('Y'); 
grid on