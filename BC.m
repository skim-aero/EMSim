function [xC,PC] = BC(xA,PA,xB,PB,PAB)
% Bar-Shalom/Campo combination 

PBA = transpose(PAB);
temp = (PA-PAB)*inv(PA+PB-PAB-PBA);

xC = xA+temp*(xB-xA);
PC = PA-temp*(PA-PBA);