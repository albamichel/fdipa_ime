[Tipo, njun, s, s0] = Modelagem3J;

r=(length(x)-1)/njun;

for i=1:njun
    q_t(:,i)=x((i-1)*r+1:i*r);
end
qf = 0;
Tf = MTH_Helic(s(njun+1,:),s0(njun+1,:),qf,'F');

[T, ~] = MTH_HelicSucess(s,s0,q_t,Tf,Tipo,njun);

Junta=1;
i=1;
for k=1:r
    J1(k,:) = [T(1,4,i) T(2,4,i)];
    i=i+3;
end

Junta=2;
i=2;
si=s(1:Junta-1,:);
s0i=s0(1:Junta-1,:);
Tf = MTH_Helic(s(Junta,:),s0(Junta,:),qf,'F');
[~, Tkk] = MTH_HelicSucess(si,s0i,q_t(1:r,1:Junta-1),Tf,Tipo,Junta-1);
for k=1:r
    J2(k,:) = [Tkk(1,4,k) Tkk(2,4,k)];
    i=i+3;
end

Junta=3;
i=3;
si=s(1:Junta-1,:);
s0i=s0(1:Junta-1,:);
Tf = MTH_Helic(s(Junta,:),s0(Junta,:),qf,'F');
[~, Tkk] = MTH_HelicSucess(si,s0i,q_t(1:r,1:Junta-1),Tf,Tipo,Junta-1);

for k=1:r
    J3(k,:) = [Tkk(1,4,k) Tkk(2,4,k)];
    i=i+3;
end

Tf = MTH_Helic(s(njun+1,:),s0(njun+1,:),qf,'F');
[T, Tk] = MTH_HelicSucess(s,s0,q_t,Tf,Tipo,njun);

%EXTRAÇÃO DAS POSIÇÕES DA FERRAMENTA (X,Y,Z)
Qx(:)=Tk(1,4,1:r);
Qy(:)=Tk(2,4,1:r);
Qz(:)=Tk(3,4,1:r);
    
k=1;
for k=1:r
clf

plot(J1(k,1), J1(k,2),'ok')
hold on
plot(J2(k,1), J2(k,2),'ok')
plot(J3(k,1), J3(k,2),'ok')
plot(Qx(1:k), Qy(1:k),'.r')
grid on

%axis('equal')   %coloca a figura no centro
axis([-10 20 -10 20])
    
line([J1(k,1), J2(k,1)],[J1(k,2), J2(k,2)])
line([J2(k,1),J3(k,1)],[J2(k,2), J3(k,2)])
line([J3(k,1), Qx(k)],[J3(k,2), Qy(k)])
pause(0.5);

end

%ARRUMAR PARA PLOTAR QUANDO R=1 (SÓ 1 PASSO)
%dq_t(1,:) = zeros(1,nvar);
delta_t = 0.1;
for i=2:r
    dq_t(i-1,:) = (q_t(i,:) - q_t(i-1,:))/delta_t;
end

figure (2) 
%posicao das juntas
subplot(2,1,1);
plot(linspace(0,r*delta_t,r),q_t);
axis([0 r*delta_t -2 2]);
xlabel('[s]');
ylabel('[rad]');
title('Posição das juntas');
legend('Junta 1', 'Junta 2','Junta 3','Location','BestOutside');

%velocidade das juntas
subplot(2,1,2);
plot(linspace(0,r*delta_t,r-delta_t),dq_t);
axis([0 r*delta_t -3 3]);
xlabel('[s]');
  ylabel('[rad/s]');
  title('Velocidade das juntas');
  legend('Junta 1', 'Junta 2','Junta 3','Location','BestOutside');

