function [q_t, T] = traptraj(qi,qf,dq_max)

dqi= [ 0 0 0 0 ];
delta_q = (qf-qi);
alpha = 0.75; % (alpha: % do tempo total com velocidade cte (dqmax) )

T = round(2*max(delta_q)/(dq_max*(1+alpha)));

T1 = (1-alpha)*T/2;
T2 = alpha*T;

for i = 1:size(qi,2)
    dqmax(i) = (delta_q(i) - dqi(i)*T1)/(T1+T2);
    %dqmax(i) = delta_q(i)/(T1);
    ddq(i) = (dqmax(i) - dqi(i)) / T1;
    %ddq(i) = 2*dqmax(i)/(T-T*alpha);
end

delta_t=0.05;
k=1;
for t = 0:delta_t:round(T)  % erro de descontinuidade na mudança de faixa (p/ alguns alpha) 
    if t <= T1
        for j=1:size(qi,2)
            q_t(k,j) = qi(j) + dqi(j)*t + 0.5*ddq(j)*t^2;
        end
        q_t1 = q_t(k,:);
    elseif t <= T1+T2 
        for j=1:size(qi,2)
            q_t(k,j) = q_t1(j) + dqmax(j)*(t-T1);
        end
        q_t2 = q_t(k,:);
    else
        for j=1:size(qi,2)
            q_t(k,j) = q_t2(j) + dqmax(j)*(t-T1-T2) - 0.5*ddq(j)*(t-T1-T2)^2;
        end
    end
    k=k+1;
end
    
end

