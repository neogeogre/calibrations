function param = compensMagnetometerSimple9(obs, init)
% compensMagnetometer does a compensation for a Magnetometer model.
% With :
% T = [t1 t2 t3
%      t4 t5 t6
%      t7 t8 t9]
%
% B = [t13; t14; t15]
%
% m = T^-1 * K * T * obs - B
%
% © Geoffrey Vincent - EPFL - 10/03/2016

% observations, all your observations goes here.
disp('==============================');
fprintf('Il y a %.0f observations \n', size(obs,1))

% Approched value, you need non zeros if the model is non linear.
t1 = init(1,1);
t2 = init(2,1);
t3 = init(3,1);
t4 = init(4,1);
t5 = init(5,1);
t6 = init(6,1);
t7 = init(7,1);
t8 = init(8,1);
t9 = init(9,1);

% Modele stochastique: matrice de covariance des observations
% Qll = eye(length(obs),length(obs));
P = eye(length(obs));
Qll = inv(P);

% Start compensation iterative
% Prealocation
it = 0;
dx = ones(size(init, 1),1);
A  = zeros(size(obs, 1), length(init));
B  = zeros(size(obs, 1), length(init)*size(obs, 1));
L  = zeros(length(obs), 1);
%tic();

while max(abs(dx))> 0.00000001 && it < 50
   
   it = it + 1;

   for i = 1:size(obs,1)
       
       % Calcul des derivees partielle
       a1 = computeDerivationt1(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a2 = computeDerivationt2(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a3 = computeDerivationt3(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a4 = computeDerivationt4(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a5 = computeDerivationt5(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a6 = computeDerivationt6(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a7 = computeDerivationt7(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a8 = computeDerivationt8(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       a9 = computeDerivationt9(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       
       mraw = computeMagnetometer(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
       l1 =  0.47 - mraw;

       % Remplissage de la matrice Jacobienne avec les derivees partielles
       A(i, 1:length(init)) = [a1, a2, a3, a4, a5, a6, a7, a8, a9];
       
%        b1 = computeDerivationm1(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
%        b2 = computeDerivationm2(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
%        b3 = computeDerivationm3(t1, t2, t3, t4, t5, t6, t7, t8, t9, obs(i,1), obs(i,2), obs(i,3));
%        
%        B(1:3) = [];
       
       % Remplissage de la matrice des ecarts apparents
       L(i, 1) = l1; % mm

   end

   % Covariance des parametres compenses Qxx
   Qxx = (A' * P * A)^-1;

   % Ajustement vector dx
   dx = Qxx * A' * P * L;

   t1 = t1 + dx(1);
   t2 = t2 + dx(2);
   t3 = t3 + dx(3);
   t4 = t4 + dx(4);
   t5 = t5 + dx(5);
   t6 = t6 + dx(6);
   t7 = t7 + dx(7);
   t8 = t8 + dx(8);
   t9 = t9 + dx(9);
   
% rank(A)
% det(A' * P * A)

end

% Correlations des parametres compenses
corrx = zeros(size(A, 2));
for i = 1:size(A, 2)
    for j = 1:size(A, 2)
        corrx(i, j) = Qxx(i,j)/sqrt(Qxx(i,i)*Qxx(j,j));
    end
end

% Create a colored plot of the correlation matrix
% figure(99)
% clf
% imagesc(abs(corrx));            
% colormap((gray)); 
% colorbar
% drawnow

V = A*dx - L; % V, matrice des residus compenses
mqo = sqrt( V'*P*V / (size(A, 1)-size(A, 2)) );

Qllcomp = A * Qxx * A';
Qvv     = Qll - Qllcomp;
sigmax  = sqrt(diag(Qxx)) * mqo;
sigmaV  = sqrt(diag(Qvv)) * mqo;
quoloc  = V ./ sigmaV; % quotient local
z = diag(Qvv*P);       % parts de redondance

% fprintf('Results after %1.0f iterations in %0.2f s : \n', it, toc());
fprintf('Results after %1.0f iterations : \n', it);
fprintf('t1 = %.7f +/- %.5f \n', t1, sigmax(1))
fprintf('t2 = %.7f +/- %.5f \n', t2, sigmax(2))
fprintf('t3 = %.7f +/- %.5f \n', t3, sigmax(3))
fprintf('t4 = %.7f +/- %.5f \n', t4, sigmax(4))
fprintf('t5 = %.7f +/- %.5f \n', t5, sigmax(5))
fprintf('t6 = %.7f +/- %.5f \n', t6, sigmax(6))
fprintf('t7 = %.7f +/- %.5f \n', t7, sigmax(7))
fprintf('t8 = %.7f +/- %.5f \n', t8, sigmax(8))
fprintf('t9 = %.7f +/- %.5f \n', t9, sigmax(9))
disp(' ');

param = [t1 t2 t3 t4 t5 t6 t7 t8 t9];
end

function m = computeMagnetometer(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

m = sqrt((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2);

end

function a1 = computeDerivationt1(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a1 = (1.0*(m1*t1-t7))*m1/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a2 = computeDerivationt2(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a2 = (1.0*(m1*t2+m2*t3-t8))*m1/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a3 = computeDerivationt3(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a3 = (1.0*(m1*t2+m2*t3-t8))*m2/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a4 = computeDerivationt4(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a4 = (1.0*(m1*t4+m2*t5+m3*t6-t9))*m1/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a5 = computeDerivationt5(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a5 = (1.0*(m1*t4+m2*t5+m3*t6-t9))*m2/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a6 = computeDerivationt6(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a6 = (1.0*(m1*t4+m2*t5+m3*t6-t9))*m3/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a7 = computeDerivationt7(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a7 = (.5*(-2*m1*t1+2*t7))/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a8 = computeDerivationt8(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a8 = (.5*(-2*m1*t2-2*m2*t3+2*t8))/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function a9 = computeDerivationt9(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

a9 = (.5*(-2*m1*t4-2*m2*t5-2*m3*t6+2*t9))/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function b1 = computeDerivationm1(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

b1 = (.5*((2*(m1*t1-t7))*t1+(2*(m1*t2+m2*t3-t8))*t2+(2*(m1*t4+m2*t5+m3*t6-t9))*t4))/ ...
     ((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function b2 = computeDerivationm2(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

b2 = (.5*((2*(m1*t2+m2*t3-t8))*t3+(2*(m1*t4+m2*t5+m3*t6-t9))*t5))/ ...
     ((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end

function b3 = computeDerivationm3(t1, t2, t3, t4, t5, t6, t7, t8, t9, m1, m2, m3)
% computeMagnetometer is the formula which is giving the squared value of the magnetic field

b3 = (1.0*(m1*t4+m2*t5+m3*t6-t9))*t6/((m1*t1-t7)^2+(m1*t2+m2*t3-t8)^2+(m1*t4+m2*t5+m3*t6-t9)^2)^.5;

end
