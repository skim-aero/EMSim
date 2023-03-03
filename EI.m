function [xC,PC] = EI(xA,PA,xB,PB)
% Ellipsoidal intersection

s = 10E-6;

% Eigenvalue decomposition of PA (First equation in eq 12a)
[EA,DA] = eig(PA);

% Second equation in eq 12a and decomposition
Temp = DA^(-0.5)*inv(EA)*PB*EA*DA^(-0.5);
[Sj,DB] = eig(Temp);

% Calculate D_gamma matrix
temp = [max(DB(1,1), 1), max(DB(2,2), 1)];
DG = diag(temp);

% Transformation matrix T and PG matrix using transformation matrix
T = EA*DA^(0.5)*Sj;
PG = T*DG*transpose(T);

% Equation 16C
H(1) = abs(DB(1,1)-1);
H(2) = abs(DB(2,2)-1);

if any(H) >= 10*s
    eta = 0;
else
    eta = s;
end

% Calculate fused mean using equation 17
xG = inv(inv(PA)+inv(PB)-2*inv(PG)+2*eta*eye(2))*...
     ((inv(PB)-inv(PG)+eta*eye(2))*xA+(inv(PA)-inv(PG)+eta*eye(2))*xB);

PC = inv(inv(PA)+inv(PB)-inv(PG));
xC = PC*(inv(PA)*xA+inv(PB)*xB-inv(PG)*xG);