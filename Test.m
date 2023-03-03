clear all;
close all;

% Mean and covariance of the estimates to be fused
PAB = [0.5 0; 0 0.5];

% xA = [1; 2];
% PA = [4 1.8; 1.8 3.5];
% 
% xB = [0.8; 1.3];
% PB = [4.5  0.5; 0.5 2.7];

% xA = [0; 0];
% PA = [4 1.8; 1.8 3.5];
% 
% xB = [0; 0];
% PB = [4.5  0.5; 0.5 2.7];

xA = [0.5; 1];
PA = [2.5 -1; -1 1.2];

xB = [2; 1];
PB = [0.8  -0.5; -0.5 4];

% xA = [0; 0];
% PA = [4 1.8; 1.8 3.5];
% 
% xB = [10; 5];
% PB = [6.5  1.5; 1.5 4.7];

% Calculate the 95% confidence line
[r_ellipseA,XA,YA] = calcEllipse(xA, PA);
[r_ellipseB,XB,YB] = calcEllipse(xB, PB);

% Bar-shalom/Campo (BC)
tic
[x_BC,P_BC] = BC(xA,PA,xB,PB,PAB);
toc
[r_ellipse_BC,X_BC,Y_BC] = calcEllipse(x_BC, P_BC);

% Convex covariance (CC)
tic
[x_CC,P_CC] = CC(xA,PA,xB,PB);
toc
[r_ellipse_CC,X_CC,Y_CC] = calcEllipse(x_CC, P_CC);

% Covariance intersection (CI) using trace minimisaiton
x(:,1) = xA;
x(:,2) = xB;

p(:,:,1) = PA;
p(:,:,2) = PB;

tic
[x_CI2,P_CI2] = fusecovint(x,p, 'trace');
toc
[r_ellipse_CI2,X_CI2,Y_CI2] = calcEllipse(x_CI2, P_CI2);

% % Coraciance intersection (CI)
% [x_CI2,P_CI2] = CI(xA,PA,xB,PB);
% [r_ellipse_CI2,X_CI2,Y_CI2] = calcEllipse(x_CI2, P_CI2);

% Set theoretic estimation (STE)
tic
[x_STE,P_STE] = STE(xA,PA,xB,PB);
toc
[r_ellipse_STE,X_STE,Y_STE] = calcEllipse(x_STE, P_STE);

% Inverse coraciance intersection (ICI)
tic
[x_ICI,P_ICI,omega] = ICI(xA,PA,xB,PB);
toc
[r_ellipse_ICI,X_ICI,Y_ICI] = calcEllipse(x_ICI, P_ICI);

% Largest ellipsoid (LE)
tic
[x_LE,P_LE] = LE(xA,PA,xB,PB);
toc
[r_ellipse_LE,X_LE,Y_LE] = calcEllipse(x_LE, P_LE);

% Internal ellipsoid approximation (IEA)
% Covariance calculation is same as LE
x0 = [2; 2]; % Starting guess
tic
[x_IEA,P_IEA] = IEA(xA,PA,xB,PB,x0);
toc
[r_ellipse_IEA,X_IEA,Y_IEA] = calcEllipse(x_IEA, P_IEA);

% Ellipsoidal intersection (EI)
tic
[x_EI,P_EI] = EI(xA,PA,xB,PB);
toc
[r_ellipse_EI,X_EI,Y_EI] = calcEllipse(x_EI, P_EI);


%% Draw the ellipsoids
plot(r_ellipseA(:,1) + XA,r_ellipseA(:,2) + YA,':','color','r', LineWidth = 2) % A
hold on
plot(r_ellipseB(:,1) + XB,r_ellipseB(:,2) + YB,':','color','b', LineWidth = 2) % B
hold on
plot(r_ellipse_CC(:,1) + X_CC,r_ellipse_CC(:,2) + Y_CC,'-','color','g', LineWidth = 2) % CC
hold on
plot(r_ellipse_BC(:,1) + X_BC,r_ellipse_BC(:,2) + Y_BC,'--','color',"#4DBEEE", LineWidth = 2) % BC
hold on
plot(r_ellipse_CI2(:,1) + X_CI2,r_ellipse_CI2(:,2) + Y_CI2,'-','color','m', LineWidth = 2) % CI
hold on
plot(r_ellipse_STE(:,1) + X_STE,r_ellipse_STE(:,2) + Y_STE,'-','color','#77AC30', LineWidth = 2) % STE
hold on
plot(r_ellipse_ICI(:,1) + X_ICI,r_ellipse_ICI(:,2) + Y_ICI,'-','color','k', LineWidth = 2) % ICI
hold on
plot(r_ellipse_LE(:,1) + X_LE,r_ellipse_LE(:,2) + Y_LE,'-','color',"#0072BD", LineWidth = 2) % LE
hold on
plot(r_ellipse_IEA(:,1) + X_IEA,r_ellipse_IEA(:,2) + Y_IEA,'-','Color',"#D95319", LineWidth = 2) % IEA
hold on
plot(r_ellipse_EI(:,1) + X_EI,r_ellipse_EI(:,2) + Y_EI,'-.','Color',"#7E2F8E", LineWidth = 2) % EI 
hold on

plot(XA, YA,'rO')
hold on
plot(XB, YB,'bO')
hold on
plot(X_CC, Y_CC,'g*')
hold on
plot(X_BC, Y_BC,'*','Color',"#4DBEEE")
hold on
plot(X_CI2, Y_CI2,'m*')
hold on
plot(X_STE, Y_STE,'*','Color',"#77AC30")
hold on
plot(X_ICI, Y_ICI,'k*')
hold on
plot(X_LE, Y_LE,'*','Color',"#0072BD")
hold on
plot(X_IEA, Y_IEA,'*','Color',"#D95319")
hold on
plot(X_EI, Y_EI,'*','Color',"#7E2F8E")
hold off

lgd = legend('A','B','CC','BC','CI','STE','ICI','LE','IEA','EI',...
       'A_{centre}','B_{centre}','CC_{centre}','BC_{centre}','CI_{centre}','STE_{centre}',...
       'ICI_{centre}','LE_{centre}','IEA_{centre}','EI_{centre}',...
       'Location','northwestoutside');
lgd.FontSize = 20;

% det_CC = det(P_CC)
% det_BC = det(P_BC)
% det_CI = det(P_CI2)
% det_STE = det(P_STE)
% det_ICI = det(P_ICI)
% det_LE = det(P_LE)
% det_IEA = det(P_IEA)
% det_EI = det(P_EI)