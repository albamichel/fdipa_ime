% %ADAPTADO PARA dq_t
clear dq_t ddq_t e_x e_y
clear Qx Qy Qz Rol Pit Yaw

ilevel = 0; % ilevel - Variables level
            % = 0 - joint position level
            % = 1 - joint velocity level
            % = 2 - joint acceleration level
            
icase = 0;  % =0 2D task, =1 3D task
            
% ROBOT PARAMETERS
%[Tipo, njun, s, s0] = Modelagem4J;
Tferr = [1 0 0 4; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%JOINT VARIABLES (r x n)
%t=linspace(0,10,r);
t=8;
%delta_t= t(2)-t(1);
r=length(t);

if ilevel == 0 % joint position level
    q_t=zeros(r,njun);
    for i=1:njun
        q_t(:,i)=x((i-1)*r+1:i*r);
    end
    
%     for i=2:r
%         dq_t(i-1,:) = (q_t(i,:) - q_t(i-1,:))/delta_t;
%     end
%         dq_t = [zeros(1,njun); dq_t];
%         
%     for i=2:r
%         ddq_t(i-1,:) = (dq_t(i,:) - dq_t(i-1,:))/delta_t;
%     end
%         ddq_t = [zeros(1,njun); ddq_t];
end

if ilevel == 1 % joint velocity level
    dq_t=zeros(r,njun);
    for i=1:njun
        dq_t(:,i)=x((i-1)*r+1:i*r);
    end
    
    q0 = [pi/12 pi/12 pi/12 pi/12];
    for i=1:njun
        q_t(:,i) = cumtrapz(t,dq_t(:,i));
        q_t(:,i) = q_t(:,i) + q0(i);
    end
end

if ilevel == 2 % joint acceleration level
    ddq_t=zeros(r,njun);
    for i=1:njun
        ddq_t(:,i)=x((i-1)*r+1:i*r);
    end
    
    dq0 = [0 0 0 0];
    for i=1:njun
        dq_t(:,i) = cumtrapz(t,ddq_t(:,i));
        dq_t(:,i) = dq_t(:,i) + dq0(i);
    end
    
    q0 = [pi/12 pi/12 pi/12 pi/12];
    for i=1:njun
        q_t(:,i) = cumtrapz(t,dq_t(:,i));
        q_t(:,i) = q_t(:,i) + q0(i);
    end
end

% JOINTS POSITION
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
    
if icase==0
    j = j+2;
end

if icase==1
    j = j+3;
end

end

[T, Tk] = MTH_HelicSucess(s,s0,q_t,Tferr,Tipo,njun);

% EXTRACT END-EFFECTOR LOCATION(px, py, pz, rx, ry, rz) FOR EACH STEP r HTM
px(:,1) = Tk(1,4,1:r);
py(:,1) = Tk(2,4,1:r);
pz(:,1) = Tk(3,4,1:r);
rx(:,1) = atan2(Tk(3,2,1:r),Tk(3,3,1:r));                            
ry(:,1) = atan2(-Tk(3,1,1:r),(Tk(3,2,1:r).^2+Tk(3,3,1:r).^2).^0.5);  
rz(:,1) = atan2(Tk(2,1,1:r),Tk(1,1,1:r));                            

% TASK EXECUTION PLOT
k=1;
for k=1:r
clf

% DESIRED PATH
pxd=rutil(:,1);
pyd=rutil(:,2);

if length(pxd) == 1
    plot(pxd,pyd,'xb');
else
    plot(pxd,pyd,'b');
end
hold on
% OBSTACLE
% Pxo = 2 - .25*cos((pi)*tt);     % posicao em X
% Pyo= 0.5  +.25*sin((pi)*tt);      % posicao em Y
% plot(Pxo,Pyo,'r:');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if floor(k/10) == k/10
% MANIPULATOR
for j=1:2:2*njun-1
    plot(J(k,j), J(k,j+1),'ok')
    hold on
end
plot(px(1:k), py(1:k),'.k')
axis([-0.5 4.5 -0.5 4.5])
for j=1:2:2*njun-3
    line([J(k,j), J(k,j+2)],[J(k,j+1), J(k,j+3)])
end
line([J(k,2*njun-1), px(k)],[J(k,2*njun), py(k)])
pause(0.5);
%end
end
   
% e_x = pdx-px;
% e_y = pdy-py;
% e_rz = rzd - rz;
% e_norm = (e_x.^2 + e_y.^2).^0.5;
    
if r > 1

    figure (2) 
    %posicao das juntas
    subplot(3,1,1);
    plot(t,q_t);
    %axis([0 r*delta_t -5 5]);
    xlabel('[s]');
    ylabel('[rad]');
    title('Posição das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');

    %velocidade das juntas
    subplot(3,1,2);
    plot(t,dq_t);
    %axis([0 r*delta_t -5 5]);
    xlabel('[s]');
    ylabel('[rad/s]');
    title('Velocidade das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
    
    %aceleração das juntas
    subplot(3,1,3);
    plot(t,ddq_t);
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
        [~,p] = csaps(t,q_t(:,k));
        fnplt(csaps(t,q_t(:,k),p/2),2,lincol(k,1:3))
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
        [~,p] = csaps(t,dq_t(:,k));
        fnplt(csaps(t,dq_t(:,k),p/2),2,lincol(k,1:3))
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
        [~,p] = csaps(t,ddq_t(:,k));
        fnplt(csaps(t,ddq_t(:,k),p/2),2,lincol(k,1:3))
        hold on
        end
    xlabel('[s]');
    ylabel('[rad/s^2]');
    title('Acelerações das juntas');
    legend('Junta 1', 'Junta 2','Junta 3','Junta 4','Location','BestOutside');
    
    figure(5)
    plot(t,normdq,'r')
    hold on
    plot(t,normddq,'b')
    xlabel('[s]');
    legend('Norma das velocidades', 'Norma das acelerações','Location','BestOutside');
end