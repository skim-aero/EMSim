function [xC,PC,omega] = ICI(xA,PA,xB,PB)
% Inverse Covariance Intersection 
%
% This function implements the ICI algorithm 
% and fuses two estimates (xA,PA) and (xB,PB).
% It provides the fusion result (xC,PC) and the 
% value of omega, which minimizes the trace of PC. 

f = @(w) trace(inv(inv(PA)+inv(PB)-inv(w*PA+(1-w)*PB)));

omega = fminbnd(f,0,1,optimset('Display','off'));
PC = inv(inv(PA)+inv(PB)-inv(omega*PA+(1-omega)*PB));

% computations of the gains
KICI = PC*(inv(PA)-omega*inv(omega*PA+(1-omega)*PB));
LICI = PC*(inv(PB)-(1-omega)*inv(omega*PA+(1-omega)*PB));

xC = KICI*xA + LICI*xB;