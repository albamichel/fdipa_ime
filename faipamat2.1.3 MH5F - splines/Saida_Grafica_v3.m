% %ADAPTADO PARA dq_t
clear dq_t ddq_t e_x e_y
clear Qx Qy Qz Rol Pit Yaw
[Tipo, njun, s, s0] = Modelagem4J;
% r = size(q_t,1);

r=(length(x)-1)/njun;
 q_t=zeros(r,njun);
% dq_t=zeros(r,njun);
% ddq_t=zeros(r,njun);
for i=1:njun
    q_t(:,i)=x((i-1)*r+1:i*r);
end

% q0 = [pi/12 pi/12 pi/12 pi/12];
% dq0 = [0 0 0 0];
% %OBTENÇÃO DE q_t INTEGRANDO dq_t
% for i=1:njun
%     dq_t(:,i) = cumtrapz(0:0.05:10,ddq_t(:,i));
%     dq_t(:,i) = dq_t(:,i) + dq0(i);
% end
% 
% for i=1:njun
%     q_t(:,i) = cumtrapz(0:0.05:10,dq_t(:,i));
%     q_t(:,i) = q_t(:,i) + q0(i);
% end

j = 1;
for Junta = 1:njun    
    if Junta==1
        i=1;
        qf = 0;
        Tf = MTH_Helic(s(njun+1,:),s0(njun+1,:),qf,'F');
        [T, ~] = MTH_HelicSucess(s,s0,q_t,Tf,Tipo,njun);
        for k=1:r
            J(k,j) = T(1,4,i);
            J(k,j+1) = T(2,4,i);
            i=i+njun;
        end
    else
        qf=0;
        Tf = MTH_Helic(s(Junta,:),s0(Junta,:),qf,'F');
        [~, Tkk] = MTH_HelicSucess(s(1:Junta-1,:),s0(1:Junta-1,:),q_t(1:r,1:Junta-1),Tf,Tipo,Junta-1);
            for k=1:r
                J(k,j) = Tkk(1,4,k);
                J(k,j+1) = Tkk(2,4,k);
            end
    end
j = j+2; %quando incluir z, alterar de 2 para 3  
end

Tf = MTH_Helic(s(njun+1,:),s0(njun+1,:),qf,'F');
[T, Tk] = MTH_HelicSucess(s,s0,q_t,Tf,Tipo,njun);

%EXTRAÇÃO DAS POSIÇÕES DA FERRAMENTA (X,Y,Z)
Qx(:) = Tk(1,4,1:r);
Qy(:) = Tk(2,4,1:r);
Qz(:) = Tk(3,4,1:r);
Yaw(:) = atan2(Tk(3,2,1:r),Tk(3,3,1:r));                            %Rx
Pit(:) = atan2(-Tk(3,1,1:r),(Tk(3,2,1:r).^2+Tk(3,3,1:r).^2).^0.5);  %Ry
Rol(:) = atan2(Tk(2,1,1:r),Tk(1,1,1:r));                            %Rz

k=1;
for k=1:r
clf

%PLOT DA TRAJETÓRIA DESEJADA
 tt=0:0.05:10;
% Px = 0.27*cos(2*pi*(sin(pi*tt/20)).^2) + 3.03905801 - 0.27;
% Py = 0.27*sin(2*pi*(sin(pi*tt/20)).^2) + 2.33195123;

Px = 2.499383626453735;
Py = 2.345266280877302;
plot(Px,Py,'xb');
hold on
%PLOT DO OBSTÁCULO
% Pxo = 2 - .25*cos((pi)*tt);     % posicao em X
% Pyo= 0.5  +.25*sin((pi)*tt);      % posicao em Y
% plot(Pxo,Pyo,'r:');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if floor(k/10) == k/10
    
for j=1:2:2*njun-1
    plot(J(k,j), J(k,j+1),'ok')
    hold on
end
plot(Qx(1:k), Qy(1:k),'.k')

%axis('equal')   %coloca a figura no centro
axis([-0.5 4.5 -0.5 4.5])

for j=1:2:2*njun-3
    line([J(k,j), J(k,j+2)],[J(k,j+1), J(k,j+3)])
end

line([J(k,2*njun-1), Qx(k)],[J(k,2*njun), Qy(k)])
pause(0.0001);

%end
end


   
% e_x = Px-Qx;
% e_y = Py-Qy;
%e_rz = (Rz - Rol).^2;

% e_norm = (e_x.^2 + e_y.^2).^0.5;
%e_ori = e_rz;
    
if r > 1
    
delta_t = 0.05;
% 
%VELOCIDADES DAS JUNTAS
%dq_t=zeros(r,njun);
for i=2:r
    dq_t(i-1,:) = (q_t(i,:) - q_t(i-1,:))/delta_t;
end
dq_t = [zeros(1,njun); dq_t];

%ACELERAÇÃO DAS JUNTAS
%ddq_t=zeros(r,njun);
for i=2:r
ddq_t(i-1,:) = (dq_t(i,:) - dq_t(i-1,:))/delta_t;
end
ddq_t = [zeros(1,njun); ddq_t];


    figure (2) 
    %posicao das juntas
    subplot(3,1,1);
    plot(tt,q_t);
    %axis([0 r*delta_t -5 5]);
    xlabel('[s]');
    ylabel('[rad]');
    title('Posição das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');

    %velocidade das juntas
    subplot(3,1,2);
    plot(tt,dq_t);
    %axis([0 r*delta_t -5 5]);
    xlabel('[s]');
    ylabel('[rad/s]');
    title('Velocidade das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
    
    %aceleração das juntas
    subplot(3,1,3);
    plot(tt,ddq_t);
    %axis([0 r*delta_t -5 5]);
    xlabel('[s]');
    ylabel('[rad/s^2]');
    title('Aceleração das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');

   
%     figure (3)
%     %erro posição normalizado
%     subplot(2,1,1);
%     plot(tt,e_norm);
%     %axis([0 r*delta_t 0 0.003]);
%     xlabel('[s]');
%     ylabel('[m]');
%     title('Erro de posição normalizado');

    %erro orientação
    % subplot(2,1,2);
    % plot(linspace(0,r*delta_t,r),e_ori);
    % xlabel('[s]');
    % ylabel('[rad]');
    % title('Erro de orientação');

lincol = ['b--'; 'g--'; 'r--'; 'k--'];

    figure(4)
    subplot(3,1,1)
        for k=1:4
        [~,p] = csaps(tt,q_t(:,k));
        fnplt(csaps(tt,q_t(:,k),p/2),2,lincol(k,1:3))
        hold on
        end
    xlabel('[s]');
    ylabel('[rad]');
    title('Posição das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
    
    
    subplot(3,1,2)
    for i=1:r
        normdq(i) = dq_t(i,:)*dq_t(i,:)';
    end
        for k=1:4
        [~,p] = csaps(tt,dq_t(:,k));
        fnplt(csaps(tt,dq_t(:,k),p/2),2,lincol(k,1:3))
        hold on
        end
    xlabel('[s]');
    ylabel('[rad/s]');
    title('Velocidades das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
   
    subplot(3,1,3)
    for i=1:r
        normddq(i) = ddq_t(i,:)*ddq_t(i,:)';
    end
        for k=1:4
        [~,p] = csaps(tt,ddq_t(:,k));
        fnplt(csaps(tt,ddq_t(:,k),p/2),2,lincol(k,1:3))
        hold on
        end
    xlabel('[s]');
    ylabel('[rad/s^2]');
    title('Acelerações das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
    
    figure(5)
    plot(tt,normdq,'r')
    hold on
    plot(tt,normddq,'b')
    xlabel('[s]');
    legend('dq^T dq', 'Norma das acelerações','Location','BestOutside');
end