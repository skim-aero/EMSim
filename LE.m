function [xC,PC] = LE(xA,PA,xB,PB)
% Largest ellipsoid method

% Calculate transformation matrix T_r and transform
[T_r,~] = eig(PA);

PA_r = T_r*PA*transpose(T_r);
PB_r = T_r*PB*transpose(T_r);

% Calculate transformation matrix T_s and transform
lambda_1 = sort(eig(PA_r));

v = [1 sqrt(lambda_1(1)/lambda_1(2))];
T_s = diag(v);

PA_sr = T_s*PA_r*transpose(T_s);
PB_sr = T_s*PB_r*transpose(T_s);

% Calculate the intersection in transformed space
[E,lambda_2] = eig(PB_sr);

k(1) = min(lambda_1(1,1), lambda_2(1,1));
k(2) = min(lambda_1(1,1), lambda_2(2,2));
D = diag(k);

E_sr = E*D*transpose(E);

% Transforming back to find fused covariance and mean
PC = inv(T_r)*inv(T_s)*E_sr*transpose(inv(T_s))*transpose(inv(T_r));
xC = inv(inv(PA)+inv(PB))*(inv(PA)*xA+inv(PB)*xB);