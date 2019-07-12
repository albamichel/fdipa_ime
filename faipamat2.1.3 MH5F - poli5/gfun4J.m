function[df,dg]=gfun4J(fun,x,indgradf,indgradg,nprob,iutil,rutil)


%
%  size(dg)= (nvar,ncstr)
%
%  indgradf=0 - df is not computed ( df=[])
%  indgradf=1 - df is computed
%
%  indgradg(i)=0 - dg(:,i) is not computed
%  indgradg(i)=1 - dg(:,i) is computed
%
%  length(indgradg)= ncstr
%
%  if indgradg(i)=0 for i=1:ncstr, take dg=zeros(nvar,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILIARY CONSTANTS
r=iutil(1);
ibeta=iutil(2);
itask=iutil(3:end);

%ROBOT PARAMETERS
[Tipo, njun, s, s0] = Modelagem4J;
Tferr = [1 0 0 4; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%JOINT VARIABLES (r x n)
q_t=zeros(r,njun);
for i=1:njun
    q_t(:,i)=x((i-1)*r+1:i*r);
end

%DESIRED TASK (pxd, pyd, pzd, rxd, ryd, rzd)
pxd=rutil(:,1);     %rxd=rutil(:,4);
pyd=rutil(:,2);     %ryd=rutil(:,5);
%pzd=rutil(:,3);    %rzd=rutil(:,6);

%COMPUTE RESULTING HTM's FOR EACH STEP r
[~, Tk] = MTH_HelicSucess(s,s0,q_t,Tferr,Tipo,njun);
    
%EXTRACT END-EFFECTOR LOCATION (px, py, pz, rx, ry, rz) FOR EACH STEP r
px(:,1) = Tk(1,4,1:r);   %rx(:,1) = atan2(Tk(3,2,1:r),Tk(3,3,1:r));
py(:,1) = Tk(2,4,1:r);   %ry(:,1) = atan2(-Tk(3,1,1:r),sqrt(Tk(3,2,1:r).^2+Tk(3,3,1:r).^2));  
%pz(:,1) = Tk(3,4,1:r);  %rz(:,1) = atan2(Tk(2,1,1:r),Tk(1,1,1:r));

%COMPUTE LOCATION ERROR (e=xd-x)
e_x = pxd-px;   %e_rx = rxd-rx;
e_y = pyd-py;   %e_ry = ryd-ry;
%e_z = pzd-pz;  %e_rz = rzd-rz;

%e_norm = sqrt(e_x.^2 + e_y.^2 + e_z.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df = []; dg=zeros(njun*r+1,0);

if indgradf==1
    df(1:length(x))=0;
    df(end)=1;
end

k=1;
j=0;
for i = 1:r
    if indgradg(i)==1
        j=j+1;
        if e_x(i) > 0 
            dg(k,j) = - (-sin(q_t(k,1)) - sin(q_t(k,1)+q_t(k,2)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(r+k,j) = - (-sin(q_t(k,1)+q_t(k,2)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(2*r+k,j) = - (-sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(3*r+k,j) = - (-sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
        else
            dg(k,j) = (-sin(q_t(k,1)) - sin(q_t(k,1)+q_t(k,2)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(r+k,j) = (-sin(q_t(k,1)+q_t(k,2)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(2*r+k,j) = (-sin(q_t(k,1)+q_t(k,2)+q_t(k,3)) - sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(3*r+k,j) = (-sin(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4)));  
        end
        k=k+1;
     end
end

k=1;
for i = r+1:2*r
    if indgradg(i)==1
        j=j+1;
        if e_y(k) > 0 
            dg(k,j) = - (cos(q_t(k,1)) + cos(q_t(k,1)+q_t(k,2)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(r+k,j) = - (cos(q_t(k,1)+q_t(k,2)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(2*r+k,j) = - (cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(3*r+k,j) = - (cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4)));
        else
            dg(k,j) = (cos(q_t(k,1)) + cos(q_t(k,1)+q_t(k,2)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(r+k,j) = (cos(q_t(k,1)+q_t(k,2)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(2*r+k,j) = (cos(q_t(k,1)+q_t(k,2)+q_t(k,3)) + cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
            dg(3*r+k,j) = (cos(q_t(k,1)+q_t(k,2)+q_t(k,3)+q_t(k,4))); 
        end
        k=k+1;
    end
end
dg(end,:) = -1;
df=df(:);

