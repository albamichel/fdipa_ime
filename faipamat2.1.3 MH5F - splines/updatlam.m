function[lambda]=updatlam(lambda0,d02)

% Updating lambda 

lambdaS=1e+5;
lamin=min(.001*d02, max(lambda0)/1000);
if lamin <= 0, lamin = .1*d02; end 
lambda=max(lambda0,lamin);
lambda=lambda(:);   
