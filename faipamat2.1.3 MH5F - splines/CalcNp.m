%Jacobiano Analítico

function [Np]= CalcNp (Caso,S0)
%syms sx sy sz
%S0 = [sx sy sz];

S=[1;0;0];
rx=[S(1); S(2); S(3); S0(2)*S(3)-S0(3)*S(2); S0(3)*S(1)-S0(1)*S(3); S0(1)*S(2)-S0(2)*S(1)];

S=[0;1;0];
ry=[S(1); S(2); S(3); S0(2)*S(3)-S0(3)*S(2); S0(3)*S(1)-S0(1)*S(3); S0(1)*S(2)-S0(2)*S(1)];

S=[0;0;1];
rz=[S(1); S(2); S(3); S0(2)*S(3)-S0(3)*S(2); S0(3)*S(1)-S0(1)*S(3); S0(1)*S(2)-S0(2)*S(1)];

S=[1;0;0];
px= [0;0;0;S(1);S(2);S(3)];

S=[0;1;0];
py= [0;0;0;S(1);S(2);S(3)];

S=[0;0;1];
pz= [0;0;0;S(1);S(2);S(3)];

if Caso=='planar'
%JA=[px(3:5) py(3:5) rz(3:5)] ;
Np=[rz(3:5) py(3:5) px(3:5)];
%JA=[rz(3:5) px(3:5) py(3:5)] ;
else
if Caso=='espacial'
Np=[rx ry rz px py pz];
end

end

