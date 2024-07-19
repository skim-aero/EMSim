function [xC,PC] = STE(xA,PA,xB,PB)
% Set theoritic estimation 

f = @(w) trace(inv(w*inv(PA)+(1-w)*inv(PB)));

omega = fminbnd(f,0,1,optimset('Display','off'));
asqre = transpose(xA-xB)*inv(inv(omega)*PA+inv(1-omega)*PB)*(xA-xB);

P0 = inv(omega*inv(PA)+(1-omega)*inv(PB));
xC = P0*(omega*inv(PA)*xA+(1-omega)*inv(PB)*xB);

PC = (1-asqre)*inv(omega*inv(PA)+(1-omega)*inv(PB));
