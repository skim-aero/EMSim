function [xC,PC] = CI(xA,PA,xB,PB)
% Covariance intersection

% f = @(w) 1/det(inv(w*inv(CA)+(1-w)*inv(CB))); 
f = @(w) trace(inv(w*inv(PA)+(1-w)*inv(PB)));

omega = fminbnd(f,0,1,optimset('Display','off'));

PC = inv(omega*inv(PA)+(1-omega)*inv(PB));
xC = PC*(omega*inv(PA)*xA+(1-omega)*inv(PB)*xB);
