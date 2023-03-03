function [xC,PC] = IEA(xA,PA,xB,PB,x0)
% Internal ellipsoid approximation

%% Calculate beta using fmincon
% f1 = @(x) transpose(x)*inv(PA)*x;
% f2 = @(x) transpose(x)*inv(PB)*x;
% 
% ceq1 = @(x) deal([], transpose(x)*inv(PB)*x-1);
% ceq2 = @(x) deal([], transpose(x)*inv(PA)*x-1);
% 
% options = optimset ('LargeScale', 'off');
% 
% [~, beta_1] = fmincon(f1, x0, [], [], [], [], [], [],...
%                       ceq1, options);
% [~, beta_2] = fmincon(f2, x0, [], [], [], [], [], [],...
%                       ceq2, options);

% Simplified beta cacluation
[x, ~] = eig(-PB*inv(PA));
beta_1 = transpose(x(:,1))*inv(PA)*x(:,1);
beta_2 = transpose(x(:,2))*inv(PB)*x(:,2);

% Calculate weight omega
omega_1 = (1-min(1,beta_2))/(1-min(1,beta_1)*min(1,beta_2));
omega_2 = (1-min(1,beta_1))/(1-min(1,beta_1)*min(1,beta_2));

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

% Calculate fused covariance and mean
xC = inv(omega_1*inv(PA)+omega_2*inv(PB))*...
        (omega_1*inv(PA)*xA+omega_2*inv(PB)*xB);

% f = @(PC) log(det(inv((1-transpose(xA)*inv(PA)*xA-transpose(xB)*inv(PB)*xB...
%            +transpose(xC)*inv(PC)*xC)*inv(omega_1*inv(PA)+omega_2*inv(PB)))));
% P0 = (PA + PB)/2;
% 
% k1 = 0.3;
% k2 = 0.7;
% c = @(PC) deal([-det([-PA*PA zeros(2) PC; zeros(2) (k1-1)*eye(2) zeros(2); PC zeros(2) -k1*eye(2)]);...
%                 -det([-PB*PB zeros(2) PC; zeros(2) (k2-1)*eye(2) zeros(2); PC zeros(2) -k2*eye(2)])], []);
% 
% options = optimset ('LargeScale', 'off');
% [PC, ~] = fmincon(f, P0, [], [], [], [], [], [],...
%                   c, options)