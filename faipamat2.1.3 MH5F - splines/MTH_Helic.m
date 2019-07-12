function T = MTH_Helic(s,s0,q,Tipo)

%syms theta d

if Tipo == 'R'
    theta=q;
    d=0;
    
a11 = cos(theta)+(s(1)^2)*(1-cos(theta)); 
a12 = s(2)*s(1)*(1-cos(theta))-s(3)*sin(theta);
a13 = s(3)*s(1)*(1-cos(theta))+s(2)*sin(theta);

a21 = s(1)*s(2)*(1-cos(theta))+s(3)*sin(theta);
a22 = cos(theta)+(s(2)^2)*(1-cos(theta));
a23 = s(3)*s(2)*(1-cos(theta))-s(1)*sin(theta);

a31 = s(1)*s(3)*(1-cos(theta))-s(2)*sin(theta);
a32 = s(2)*s(3)*(1-cos(theta))+s(1)*sin(theta);
a33 = cos(theta)+(s(3)^2)*(1-cos(theta));

qx = d*s(1)+s0(1)*(1-a11)-s0(2)*a12-s0(3)*a13;
qy = d*s(2)-s0(1)*a21+s0(2)*(1-a22)-s0(3)*a23;
qz = d*s(3)-s0(1)*a31-s0(2)*a32+s0(3)*(1-a33);

T = [a11 a12 a13 qx;
     a21 a22 a23 qy;
     a31 a32 a33 qz;
       0   0   0  1];
    
    %T = eval(T);
    
elseif Tipo == 'P'
    theta=0;
    d=q;
    
a11 = cos(theta)+(s(1)^2)*(1-cos(theta)); 
a12 = s(2)*s(1)*(1-cos(theta))-s(3)*sin(theta);
a13 = s(3)*s(1)*(1-cos(theta))+s(2)*sin(theta);

a21 = s(1)*s(2)*(1-cos(theta))+s(3)*sin(theta);
a22 = cos(theta)+(s(2)^2)*(1-cos(theta));
a23 = s(3)*s(2)*(1-cos(theta))-s(1)*sin(theta);

a31 = s(1)*s(3)*(1-cos(theta))-s(2)*sin(theta);
a32 = s(2)*s(3)*(1-cos(theta))+s(1)*sin(theta);
a33 = cos(theta)+(s(3)^2)*(1-cos(theta));

qx = d*s(1)+s0(1)*(1-a11)-s0(2)*a12-s0(3)*a13;
qy = d*s(2)-s0(1)*a21+s0(2)*(1-a22)-s0(3)*a23;
qz = d*s(3)-s0(1)*a31-s0(2)*a32+s0(3)*(1-a33);

T = [a11 a12 a13 qx;
     a21 a22 a23 qy;
     a31 a32 a33 qz;
       0   0   0  1];
   
    %T = eval(T);

elseif Tipo == 'F'
    theta=q;
    d=0;
    
a11 = cos(theta)+(s(1)^2)*(1-cos(theta)); 
a12 = s(2)*s(1)*(1-cos(theta))-s(3)*sin(theta);
a13 = s(3)*s(1)*(1-cos(theta))+s(2)*sin(theta);

a21 = s(1)*s(2)*(1-cos(theta))+s(3)*sin(theta);
a22 = cos(theta)+(s(2)^2)*(1-cos(theta));
a23 = s(3)*s(2)*(1-cos(theta))-s(1)*sin(theta);

a31 = s(1)*s(3)*(1-cos(theta))-s(2)*sin(theta);
a32 = s(2)*s(3)*(1-cos(theta))+s(1)*sin(theta);
a33 = cos(theta)+(s(3)^2)*(1-cos(theta));

qx = s0(1);
qy = s0(2);
qz = s0(3);

T = [a11 a12 a13 qx;
     a21 a22 a23 qy;
     a31 a32 a33 qz;
       0   0   0  1];
    
    %T(1:3,4) = [qx; qy; qz]';
    %T = eval(T);
end

end

