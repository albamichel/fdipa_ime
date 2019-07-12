%Cálculo de Ns

function [Ns,MatAux]= CalcNs (T,s0ref)

MatAux(:,1) = [s0ref(1,:)'; 1];
Ns(:,1) = [1 MatAux(2,1) -MatAux(1,1)];

MatAux(:,2) = T(:,:,1)*[s0ref(2,:)'; 1];
Ns(:,2) = [1 MatAux(2,2) -MatAux(1,2)];

MatAux(:,3) = T(:,:,1)*T(:,:,2)*[s0ref(3,:)'; 1];
Ns(:,3) = [1 MatAux(2,3) -MatAux(1,3)];

MatAux(:,4) = T(:,:,1)*T(:,:,2)*T(:,:,3)*[s0ref(4,:)'; 1];
Ns(:,4) = [1 MatAux(2,4) -MatAux(1,4)];

MatAux(:,5) = T(:,:,1)*T(:,:,2)*T(:,:,3)*T(:,:,4)*[s0ref(5,:)'; 1];
Ns(:,5) = [1 MatAux(2,5) -MatAux(1,5)];

MatAux(:,6) = T(:,:,1)*T(:,:,2)*T(:,:,3)*T(:,:,4)*T(:,:,5)*[s0ref(6,:)'; 1];
Ns(:,6) = [1 MatAux(2,6) -MatAux(1,6)];

end