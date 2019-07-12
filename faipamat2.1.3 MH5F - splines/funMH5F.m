
function[f,g]=fun4J(x,indf,indg,nprob,iutil,rutil)

%  fun32 - computes f and g 
% 
%  f is the objective function
%  g(ncstr) is a vector with the response contraints ordered as follows:
% 
%  g(1),...,g(neqlin)     - Linear equality constraints
%  g(neqlin+1),...,g(neq) - Nonlinear equality constraints
%  g(neq+1),...,g(ncstr)  - Inequality constraints (it does not
%                           include box constraints)
%
%  indf=0 - f is not computed ( Take f=[])
%  indf=1 - f is computed
%
%  indg(i)=0 - g(i) is not computed
%  indg(i)=1 - g(i) is computed
%
%  length(indg)= ncstr
%
%  if indg(i)=0 for i=1:ncstr, take g=[]
%
%
%  
%  iutil - Integer utility vector, employed to store data in
%                  probXX, funXX and gfunXX.
% 
%  rutil - Real utility vector, employed to 
%                  store data in probXX, funXX and gfunXX.
%
     
f=[];g=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILIARY CONSTANTS
r=iutil(1);
ibeta=iutil(2);
itask=iutil(3:end);
t=linspace(0,2,r);

%ROBOT PARAMETERS
[Tipo, njun, s, s0] = ModelagemMH5F;
Tferr = [0 0 1 .585; 0 -1 0 0; 1 0 0 .350; 0 0 0 1];

%JOINT VARIABLES (r x n)
% q_t=zeros(r,njun);
% for i=1:njun
%     q_t(:,i)=x((i-1)*r+1:i*r);
% end

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

%t=0:.01:1;
%r=length(t);

q1(:,1) = a(1,1) + a(1,2)*t + a(1,3)*t.^2 + a(1,4)*t.^3 + a(1,5)*t.^4 + a(1,6)*t.^5;
q2(:,1) = a(2,1) + a(2,2)*t + a(2,3)*t.^2 + a(2,4)*t.^3 + a(2,5)*t.^4 + a(2,6)*t.^5;
q3(:,1) = a(3,1) + a(3,2)*t + a(3,3)*t.^2 + a(3,4)*t.^3 + a(3,5)*t.^4 + a(3,6)*t.^5;
q4(:,1) = a(4,1) + a(4,2)*t + a(4,3)*t.^2 + a(4,4)*t.^3 + a(4,5)*t.^4 + a(4,6)*t.^5;
q5(:,1) = a(5,1) + a(5,2)*t + a(5,3)*t.^2 + a(5,4)*t.^3 + a(5,5)*t.^4 + a(5,6)*t.^5;
q6(:,1) = a(6,1) + a(6,2)*t + a(6,3)*t.^2 + a(6,4)*t.^3 + a(6,5)*t.^4 + a(6,6)*t.^5;

q_t = [q1 q2 q3 q4 q5 q6];

%DESIRED TASK (pxd, pyd, pzd, rxd, ryd, rzd)
pxd=rutil(:,1);
pyd=rutil(:,2);
pzd=rutil(:,3);

%COMPUTE RESULTING HTM's FOR EACH STEP r
[~, Tk] = MTH_HelicSucess(s,s0,q_t,Tferr,Tipo,njun);
    
%EXTRACT END-EFFECTOR LOCATION (px, py, pz, rx, ry, rz) FOR EACH STEP r
px(:,1) = Tk(1,4,1:r);   %rx(:,1) = atan2(Tk(3,2,1:r),Tk(3,3,1:r));
py(:,1) = Tk(2,4,1:r);   %ry(:,1) = atan2(-Tk(3,1,1:r),sqrt(Tk(3,2,1:r).^2+Tk(3,3,1:r).^2));  
pz(:,1) = Tk(3,4,1:r);   %rz(:,1) = atan2(Tk(2,1,1:r),Tk(1,1,1:r));

%COMPUTE LOCATION ERROR (e=xd-x)
e_x = pxd-px;   %e_rx = rxd-rx;
e_y = pyd-py;   %e_ry = ryd-ry;
e_z = pzd-pz;   %e_rz = rzd-rz;

%e_norm = sqrt(e_x.^2 + e_y.^2 + e_z.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT POSITIONS
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

%     if icase==0
%        j = j+2;
%     end

%     if icase==1
         j = j+3;
%     end

end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:r
    delta_x(i-1,1) = sqrt((px(i)-px(i-1))^2 + (py(i)-py(i-1))^2 + (pz(i)-pz(i-1))^2);
end

if indf==1
    f = sum(delta_x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dq1(:,1) = a(1,2) + 2*a(1,3)*t + 3*a(1,4)*t.^2 + 4*a(1,5)*t.^3 + 5*a(1,6)*t.^4;
dq2(:,1) = a(2,2) + 2*a(2,3)*t + 3*a(2,4)*t.^2 + 4*a(2,5)*t.^3 + 5*a(2,6)*t.^4;
dq3(:,1) = a(3,2) + 2*a(3,3)*t + 3*a(3,4)*t.^2 + 4*a(3,5)*t.^3 + 5*a(3,6)*t.^4;
dq4(:,1) = a(4,2) + 2*a(4,3)*t + 3*a(4,4)*t.^2 + 4*a(4,5)*t.^3 + 5*a(4,6)*t.^4;
dq5(:,1) = a(5,2) + 2*a(5,3)*t + 3*a(5,4)*t.^2 + 4*a(5,5)*t.^3 + 5*a(5,6)*t.^4;
dq6(:,1) = a(6,2) + 2*a(6,3)*t + 3*a(6,4)*t.^2 + 4*a(6,5)*t.^3 + 5*a(6,6)*t.^4;

ddq1(:,1) = 2*a(1,3) + 6*a(1,4)*t + 12*a(1,5)*t.^2 + 20*a(1,6)*t.^3;
ddq2(:,1) = 2*a(2,3) + 6*a(2,4)*t + 12*a(2,5)*t.^2 + 20*a(2,6)*t.^3;
ddq3(:,1) = 2*a(3,3) + 6*a(3,4)*t + 12*a(3,5)*t.^2 + 20*a(3,6)*t.^3;
ddq4(:,1) = 2*a(4,3) + 6*a(4,4)*t + 12*a(4,5)*t.^2 + 20*a(4,6)*t.^3;
ddq5(:,1) = 2*a(5,3) + 6*a(5,4)*t + 12*a(5,5)*t.^2 + 20*a(5,6)*t.^3;
ddq6(:,1) = 2*a(6,3) + 6*a(6,4)*t + 12*a(6,5)*t.^2 + 20*a(6,6)*t.^3;

% j=0;
% if indg(1)==1
%     j=j+1;
%     g(j) = abs(px(1) - 1.5) - 0.001;
% end
% 
% if indg(2)==1
%     j=j+1;
%     g(j) = abs(py(1) - 1.5) - 0.001;
% end

j=0;
if indg(1)==1
    j=j+1;
    g(j) = abs(q1(1)) - 0.0001;
end

if indg(2)==1
    j=j+1;
    g(j) = abs(q2(1)) - 0.0001;
end

if indg(3)==1
    j=j+1;
    g(j) = abs(q3(1)) - 0.0001;
end

if indg(4)==1
    j=j+1;
    g(j) = abs(q4(1)) - 0.0001;
end

if indg(5)==1
    j=j+1;
    g(j) = abs(q5(1)) - 0.0001;
end
if indg(6)==1
    j=j+1;
    g(j) = abs(q6(1)) - 0.0001;
end

if indg(7)==1
    j=j+1;
    g(j) = abs(px(end) - pxd) - 0.001;
end

if indg(8)==1
    j=j+1;
    g(j) = abs(py(end) - pyd) - 0.001;
end

if indg(9)==1
    j=j+1;
    g(j) = abs(pz(end) - pzd) - 0.001;
end

if indg(10)==1
    j=j+1;
    g(j) = abs(dq1(1)) - 0.001;
end

if indg(11)==1
    j=j+1;
    g(j) = abs(dq2(1)) - 0.001;
end

if indg(12)==1
    j=j+1;
    g(j) = abs(dq3(1)) - 0.001;
end

if indg(13)==1
    j=j+1;
    g(j) = abs(dq4(1)) - 0.001;
end

if indg(14)==1
    j=j+1;
    g(j) = abs(dq5(1)) - 0.001;
end
if indg(15)==1
    j=j+1;
    g(j) = abs(dq6(1)) - 0.001;
end

if indg(16)==1
    j=j+1;
    g(j) = abs(dq1(end)) - 0.001;
end

if indg(17)==1
    j=j+1;
    g(j) = abs(dq2(end)) - 0.001;
end

if indg(18)==1
    j=j+1;
    g(j) = abs(dq3(end)) - 0.001;
end

if indg(19)==1
    j=j+1;
    g(j) = abs(dq4(end)) - 0.001;
end

if indg(20)==1
    j=j+1;
    g(j) = abs(dq5(end)) - 0.001;
end
if indg(21)==1
    j=j+1;
    g(j) = abs(dq6(end)) - 0.001;
end

if indg(22)==1
    j=j+1;
    g(j) = abs(ddq1(1)) - 0.001;
end

if indg(23)==1
    j=j+1;
    g(j) = abs(ddq2(1)) - 0.001;
end

if indg(24)==1
    j=j+1;
    g(j) = abs(ddq3(1)) - 0.001;
end

if indg(25)==1
    j=j+1;
    g(j) = abs(ddq4(1)) - 0.001;
end

if indg(26)==1
    j=j+1;
    g(j) = abs(ddq5(1)) - 0.001;
end
if indg(27)==1
    j=j+1;
    g(j) = abs(ddq6(1)) - 0.001;
end

if indg(28)==1
    j=j+1;
    g(j) = abs(ddq1(end)) - 0.001;
end

if indg(29)==1
    j=j+1;
    g(j) = abs(ddq2(end)) - 0.001;
end

if indg(30)==1
    j=j+1;
    g(j) = abs(ddq3(end)) - 0.001;
end

if indg(31)==1
    j=j+1;
    g(j) = abs(ddq4(end)) - 0.001;
end

if indg(32)==1
    j=j+1;
    g(j) = abs(ddq5(end)) - 0.001;
end
if indg(33)==1
    j=j+1;
    g(j) = abs(ddq6(end)) - 0.001;
end

% if indg(13)==1
%     j=j+1;
%     g(j) = abs(ddq1(1)) - 0.001;
% end
% 
% if indg(14)==1
%     j=j+1;
%     g(j) = abs(ddq2(1)) - 0.001;
% end
% 
% if indg(15)==1
%     j=j+1;
%     g(j) = abs(ddq3(1)) - 0.001;
% end
% 
% if indg(16)==1
%     j=j+1;
%     g(j) = abs(ddq4(1)) - 0.001;
% end
% 
% if indg(17)==1
%     j=j+1;
%     g(j) = abs(ddq1(end)) - 0.001;
% end
% 
% if indg(18)==1
%     j=j+1;
%     g(j) = abs(ddq2(end)) - 0.001;
% end
% 
% if indg(19)==1
%     j=j+1;
%     g(j) = abs(ddq3(end)) - 0.001;
% end
% 
% if indg(20)==1
%     j=j+1;
%     g(j) = abs(ddq4(end)) - 0.001;
% end

obs=[ .585; 0; .2 ];
k=1;
for i=34:33+r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs-[px(k); py(k); pz(k)]);
    end
    k=k+1;
end

k=1;
for i=34+r:33+2*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,1:3));
    end
    k=k+1;
end

k=1;
for i=34+2*r:33+3*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,4:6));
    end
    k=k+1;
end

k=1;
for i=34+3*r:33+4*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,7:9));
    end
    k=k+1;
end

k=1;
for i=34+4*r:33+5*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,10:12));
    end
    k=k+1;
end

k=1;
for i=34+5*r:33+6*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,13:15));
    end
    k=k+1;
end

k=1;
for i=34+6*r:33+7*r
    if indg(i)==1
        j=j+1; 
        g(j) = 0.1 - norm(obs'-p(k,16:18));
    end
    k=k+1;
end

k=1;
for i=34+7*r:33+8*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q1(k)) - pi;
    end
    k=k+1;
end

k=1;
for i=34+8*r:33+9*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q2(k)) - pi;
    end
    k=k+1;
end

k=1;
for i=34+9*r:33+10*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q3(k)) - pi;
    end
    k=k+1;
end

k=1;
for i=34+10*r:33+11*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q4(k)) - pi;
    end
    k=k+1;
end

k=1;
for i=34+11*r:33+12*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q5(k)) - pi;
    end
    k=k+1;
end

k=1;
for i=34+12*r:33+13*r
    if indg(i)==1
        j=j+1; 
        g(j) = norm(q6(k)) - pi;
    end
    k=k+1;
end
g=g(:);