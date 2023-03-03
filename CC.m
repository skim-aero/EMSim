function [xC,PC] = CC(xA,PA,xB,PB)
% Convex comnbination

PC = inv(inv(PA)+inv(PB));
xC = PC*(inv(PA)*xA+inv(PB)*xB);