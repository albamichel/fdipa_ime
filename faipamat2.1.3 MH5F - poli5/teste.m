clear all
qi = [ 0 0 0 0];

qf = [0.070470864764514
      0.380255467943407
      0.548562236313176
      0.510636255508920]';
  
dqmax = 0.075;
  [q_t,Tt] = traptraj(qi,qf,dqmax);
  
  xt=q_t(:);