% %ADAPTADO PARA dq_t
clear dq_t ddq_t e_x e_y px py pz rx ry rz
clear Qx Qy Qz Rol Pit Yaw

ilevel = 0; % ilevel - Variables level
            % = 0 - joint position level
            % = 1 - joint velocity level
            % = 2 - joint acceleration level
            
icase = 1;  % =0 2D task, =1 3D task
            
% ROBOT PARAMETERS
[Tipo, njun, s, s0] = ModelagemMH5F;
Tferr = [0 0 1 .585; 0 -1 0 0; 1 0 0 .350; 0 0 0 1];

%JOINT VARIABLES (r x n)
%t=linspace(0,1,r);
delta_t= t(2)-t(1);
%t=1;
a(1,1) = x(1); %junta1
a(1,2) = x(2);
a(1,3) = x(3);
a(1,4) = x(4);
a(1,5) = x(5);
a(1,6) = x(6);

a(2,1) = x(7); %junta2
a(2,2) = x(8);
a(2,3) = x(9);
a(2,4) = x(10);
a(2,5) = x(11);
a(2,6) = x(12);

a(3,1) = x(13); %junta3
a(3,2) = x(14);
a(3,3) = x(15);
a(3,4) = x(16);
a(3,5) = x(17);
a(3,6) = x(18);

a(4,1) = x(19); %junta4
a(4,2) = x(20);
a(4,3) = x(21);
a(4,4) = x(22);
a(4,5) = x(23);
a(4,6) = x(24);

a(5,1) = x(25); %junta5
a(5,2) = x(26);
a(5,3) = x(27);
a(5,4) = x(28);
a(5,5) = x(29);
a(5,6) = x(30);

a(6,1) = x(31); %junta6
a(6,2) = x(32);
a(6,3) = x(33);
a(6,4) = x(34);
a(6,5) = x(35);
a(6,6) = x(36);

%t=0:.1:1;
%r=length(t);

q1(:,1) = a(1,1) + a(1,2)*t + a(1,3)*t.^2 + a(1,4)*t.^3 + a(1,5)*t.^4 + a(1,6)*t.^5;
q2(:,1) = a(2,1) + a(2,2)*t + a(2,3)*t.^2 + a(2,4)*t.^3 + a(2,5)*t.^4 + a(2,6)*t.^5;
q3(:,1) = a(3,1) + a(3,2)*t + a(3,3)*t.^2 + a(3,4)*t.^3 + a(3,5)*t.^4 + a(3,6)*t.^5;
q4(:,1) = a(4,1) + a(4,2)*t + a(4,3)*t.^2 + a(4,4)*t.^3 + a(4,5)*t.^4 + a(4,6)*t.^5;
q5(:,1) = a(5,1) + a(5,2)*t + a(5,3)*t.^2 + a(5,4)*t.^3 + a(5,5)*t.^4 + a(5,6)*t.^5;
q6(:,1) = a(6,1) + a(6,2)*t + a(6,3)*t.^2 + a(6,4)*t.^3 + a(6,5)*t.^4 + a(6,6)*t.^5;

q_t=[q1 q2 q3 q4 q5 q6];

dq1(:,1) = a(1,2) + 2*a(1,3)*t + 3*a(1,4)*t.^2 + 4*a(1,5)*t.^3 + 5*a(1,6)*t.^4;
dq2(:,1) = a(2,2) + 2*a(2,3)*t + 3*a(2,4)*t.^2 + 4*a(2,5)*t.^3 + 5*a(2,6)*t.^4;
dq3(:,1) = a(3,2) + 2*a(3,3)*t + 3*a(3,4)*t.^2 + 4*a(3,5)*t.^3 + 5*a(3,6)*t.^4;
dq4(:,1) = a(4,2) + 2*a(4,3)*t + 3*a(4,4)*t.^2 + 4*a(4,5)*t.^3 + 5*a(4,6)*t.^4;
dq5(:,1) = a(5,2) + 2*a(5,3)*t + 3*a(5,4)*t.^2 + 4*a(5,5)*t.^3 + 5*a(5,6)*t.^4;
dq6(:,1) = a(6,2) + 2*a(6,3)*t + 3*a(6,4)*t.^2 + 4*a(6,5)*t.^3 + 5*a(6,6)*t.^4;

dq_t=[dq1 dq2 dq3 dq4 dq5 dq6];

ddq1(:,1) = 2*a(1,3) + 6*a(1,4)*t + 12*a(1,5)*t.^2 + 20*a(1,6)*t.^3;
ddq2(:,1) = 2*a(2,3) + 6*a(2,4)*t + 12*a(2,5)*t.^2 + 20*a(2,6)*t.^3;
ddq3(:,1) = 2*a(3,3) + 6*a(3,4)*t + 12*a(3,5)*t.^2 + 20*a(3,6)*t.^3;
ddq4(:,1) = 2*a(4,3) + 6*a(4,4)*t + 12*a(4,5)*t.^2 + 20*a(4,6)*t.^3;
ddq5(:,1) = 2*a(5,3) + 6*a(5,4)*t + 12*a(5,5)*t.^2 + 20*a(5,6)*t.^3;
ddq6(:,1) = 2*a(6,3) + 6*a(6,4)*t + 12*a(6,5)*t.^2 + 20*a(6,6)*t.^3;

ddq_t=[ddq1 ddq2 ddq3 ddq4 ddq5 ddq6];

dddq1(:,1) = 6*a(1,4) + 24*a(1,5)*t + 60*a(1,6)*t.^2;
dddq2(:,1) = 6*a(2,4) + 24*a(2,5)*t + 60*a(2,6)*t.^2;
dddq3(:,1) = 6*a(3,4) + 24*a(3,5)*t + 60*a(3,6)*t.^2;
dddq4(:,1) = 6*a(4,4) + 24*a(4,5)*t + 60*a(4,6)*t.^2;
dddq5(:,1) = 6*a(5,4) + 24*a(5,5)*t + 60*a(5,6)*t.^2;
dddq6(:,1) = 6*a(6,4) + 24*a(6,5)*t + 60*a(6,6)*t.^2;

dddq_t=[dddq1 dddq2 dddq3 dddq4 dddq5 dddq6];

if ilevel == 0 % joint position level
%     q_t=zeros(r,njun);
%     for i=1:njun
%         q_t(:,i)=x((i-1)*r+1:i*r);
%     end

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
            J(k,j+2) = T(3,4,i);
            i=i+njun;
        end
    else
        qf=0;
        Tf = MTH_Helic(s(Junta,:),s0(Junta,:),qf,'F');
        [~, Tkk] = MTH_HelicSucess(s(1:Junta-1,:),s0(1:Junta-1,:),q_t(1:r,1:Junta-1),Tf,Tipo,Junta-1);
            for k=1:r
                J(k,j) = Tkk(1,4,k);
                J(k,j+1) = Tkk(2,4,k);
                J(k,j+2) = Tkk(3,4,k);                
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

%CENTRO DOS ELOS
p1 = J(:,1:3);
p2 = J(:,4:6);
p3 = J(:,7:9);
p4 = J(:,10:12);
p5 = J(:,13:15);
p6 = J(:,16:18);
p = [p1 p2 p3 p4 p5 p6];

p1a = (p1+p2)/2;
p2a = (p2+p3)/2;
p3a = (p3+p4)/2;
p4a = (p4+p5)/2;
p5a = (p5+p6)/2;
p6a = [p6(:,1)+px(:) p6(:,2)+py(:) p6(:,3)+pz(:)]/2;
pa = [p1a p2a p3a p4a p5a p6a];

% TASK EXECUTION PLOT
k=1;
figure(1)

figcount=1;
xlegcount=1;
x_leg = ['0.0s'; '0.3s'; '4.0s'; '0.5s'; '0.6s'; '1.0s'];

pxd=rutil(:,1);
pyd=rutil(:,2);
pzd=rutil(:,3);

for k=1:r
clf

% DESIRED PATH
if length(pxd) == 1
    plot3(pxd,pyd,pzd,'xb');
else
    plot3(pxd,pyd,pzd,'r');
end
hold on
axis([0 .8 -.4 .4 -.1 .7])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if (floor(k/10) == k/10) || k==1
%if k==1 || k==31 || k==41 || k==51 || k==61  || k==101
%subplot(2,3,figcount)

% OBSTACLE
%  tt=0:0.0001:2;
% 
%  Pxo = 1.75 - .3*cos((pi)*tt);     % posicao em X
%  Pyo = 1.75 + .3*sin((pi)*tt);      % posicao em Y
%  plot(Pxo,Pyo,'r');
    obs=[ .585; 0; .2 ];
    RR=0.1;
    [aa,bb,cc]=sphere;
    surf(aa*RR+obs(1), bb*RR+obs(2), cc*RR+obs(3))
    hold on
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% MANIPULADOR
for j=1:3:3*njun-2
    plot3(J(k,j), J(k,j+1), J(k,j+2),'ok')
    hold on
    grid on
end
plot3(px(1:k), py(1:k), pz(1:k),'.b')

plot3(p1a(k,1),p1a(k,2),p1a(k,3),'hr')
plot3(p2a(k,1),p2a(k,2),p2a(k,3),'hr')
plot3(p3a(k,1),p3a(k,2),p3a(k,3),'hr')
plot3(p4a(k,1),p4a(k,2),p4a(k,3),'hr')
plot3(p5a(k,1),p5a(k,2),p5a(k,3),'hr')
plot3(p6a(k,1),p6a(k,2),p6a(k,3),'hr')

%axis([-1 3 -.5 3])
%xlabel(x_leg(xlegcount,:),'fontsize',12);
%xlegcount=xlegcount+1;
%ylabel('y [m]');

for j=1:3:3*njun-4
    %line([J(k,j), J(k,j+2)],[J(k,j+1), J(k,j+3)],'LineWidth',2,'Color',[(1/size(x_iter,1))^.95*kcolor 0 0])
        line([J(k,j), J(k,j+3)],[J(k,j+1), J(k,j+4)],[J(k,j+2), J(k,j+5)],'LineWidth',2,'Color',[0 0 0])
end
%line([J(k,2*njun-1), px(k)],[J(k,2*njun), py(k)],'LineWidth',2,'Color',[(1/size(x_iter,1))^.95*kcolor 0 0])
line([J(k,3*njun-2), px(k)],[J(k,3*njun-1), py(k)],[J(k,3*njun), pz(k)],'LineWidth',2,'Color',[0 0 0])
pause(0.01);


% if k==1
%     line([0, px(k)], [0, py(k)]);
%     hold on
% else
%    line([px(k), px(k-1)], [py(k), py(k-1)])
%    hold on
% end
% hold on
figcount=figcount+1;
%end

end
   
e_x = pxd-px;
e_y = pyd-py;
e_z = pzd-pz;
erro=[e_x e_y e_z];

%e_rz = rzd - rz;
e_norm = sqrt(e_x.^2 + e_y.^2 + e_z.^2);

lincol = ['b: '; 'r-.'; 'y- ';'k--';'k--';'k--'];
lincol2 = ['k- '; 'r- '; 'g- ';'b- ';'y- ';'m- '];
% 
% subplot(2,2,2)
%     for k=1:2
%     [~,p] = csaps(t,erro(:,k));
%     fnplt(csaps(t,erro(:,k),p/2),2,'LineWidth',1,lincol(k,1:3))
%     hold on
%     end
%     
% % for i=1:2
% %     plot(tt,erro(:,i),lincol2(i,1:3),'LineWidth',0.5);
% %     hold on
% % end
% legend('e_x','e_y','Location','south','Orientation','horizontal');
% xlabel({'[s]','(b)'});
% ylabel('[m]');
% 
% subplot(2,2,3)
%     for k=1:njun
%     [~,p] = csaps(t,q_t(:,k));
%     fnplt(csaps(t,q_t(:,k),p/2),2,'LineWidth',1,lincol(k,1:3))
%     hold on
%     end
%     xlabel('(c)');
%     %ylabel('[rad]');
%     %title('Joint Position');
%     legend('\theta_1', '\theta_2','\theta_3','Location','north','Orientation','horizontal');
% 
% subplot(2,2,4)
%     for i=1:r
%         normdq(i) = dq_t(i,:)*dq_t(i,:)';
%     end
%         for k=1:njun
%         [~,p] = csaps(t,dq_t(:,k));
%         fnplt(csaps(t,dq_t(:,k),p/2),2,'LineWidth',1,lincol(k,1:3))
%         hold on
%         end
%     xlabel('(d)');
%     %ylabel('[rad/s]');
%     %title('Joint Velocity');
%     legend({'$\dot\theta_1$', '$\dot\theta_2$','$\dot\theta_3$'},'Location','north','Orientation','horizontal','Interpreter','latex');
%     
if r > 1

    figure (2) 
    %posicao das juntas
    %subplot(2,1,1);
    for k=1:njun
        plot(t,q_t(:,k),lincol2(k,1:3),'LineWidth',1)
        hold on
        grid on
    end
    %axis([0 t(end) -2 3.5]);
    xlabel('[s]','fontsize',12);
    ylabel('[rad]','fontsize',12);
    %title('Posição das juntas','fontsize',12);
    legend({'\theta_1', '\theta_2','\theta_3','\theta_4','\theta_5','\theta_6'},'fontsize',12,'Location','NorthEast','Orientation','horizontal');
    
    figure (3)
    %velocidade das juntas
    %subplot(2,1,2);
    for k=1:njun
        plot(t,dq_t(:,k),lincol2(k,1:3),'LineWidth',1)
        hold on
        grid on
    end
    %axis([0 t(end) -2.5 2.5]);
    xlabel('[s]','fontsize',12);
    ylabel('[rad/s]','fontsize',12);
    %title('Velocidade das juntas','fontsize',12);
    legend({'$\dot\theta_1$', '$\dot\theta_2$','$\dot\theta_3$','$\dot\theta_4$','$\dot\theta_5$','$\dot\theta_6$'},'fontsize',12,'Location','NorthEast','Orientation','horizontal','Interpreter','latex');
    
    figure (4)
    %aceleracao das juntas
    %subplot(2,1,2);
    for k=1:njun
        plot(t,ddq_t(:,k),lincol2(k,1:3),'LineWidth',1)
        hold on
        grid on
    end
    %axis([0 t(end) -2.5 2.5]);
    xlabel('[s]','fontsize',12);
    ylabel('[rad/s2]','fontsize',12);
    %title('Velocidade das juntas','fontsize',12);
    legend({'$\ddot\theta_1$', '$\ddot\theta_2$','$\ddot\theta_3$','$\ddot\theta_4$','$\ddot\theta_5$','$\ddot\theta_6$'},'fontsize',12,'Location','NorthEast','Orientation','horizontal','Interpreter','latex');
    
 figure (5)
    %jerk das juntas
    %subplot(2,1,2);
    for k=1:njun
        plot(t,dddq_t(:,k),lincol2(k,1:3),'LineWidth',1)
        hold on
        grid on
    end
    %axis([0 t(end) -2.5 2.5]);
    xlabel('[s3]','fontsize',12);
    ylabel('[rad/s2]','fontsize',12);
    %title('Velocidade das juntas','fontsize',12);
    legend({'$\stackrel{...}{\theta_1}$', '$\stackrel{...}{\theta_2}$','$\stackrel{...}{\theta_3}$','$\stackrel{...}{\theta_4}$','$\stackrel{...}{\theta_5}$','$\stackrel{...}{\theta_6}$'},'fontsize',12,'Location','NorthEast','Orientation','horizontal','Interpreter','latex');
   
   
%     figure (3)
%     %erro posição normalizado
%     subplot(2,1,1);
%     plot(t,e_norm,'LineWidth',1);
%     hold on
%     %axis([0 r*delta_t 0 0.003]);
%     %xlabel('[s]','fontsize',12);
%     ylabel('[m]','fontsize',12);
%     title('Erro de posição normalizado','fontsize',12);
% 
%     %erro orientação
%     subplot(2,1,2);
%     plot(t,e_rz,'LineWidth',1);
%     hold on
%     %xlabel('[s]','fontsize',12);
%     ylabel('[rad]','fontsize',12);
%     title('Erro de orientação','fontsize',12);

% NORMA DAS VELOCIDADES E ACELERAÇÕES ***********************************
% figure(5)
%     plot(t,normdq,'r')
%     hold on
%     plot(t,normddq,'b')
%     xlabel('[s]');
%     legend('Norma das velocidades', 'Norma das acelerações','Location','BestOutside');
end