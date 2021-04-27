function param = compens6param(obs, init)
% compens9param does a leat mean square for computing bx, bxy,
% bxz, byz, in respect of fx^2 + fy^2 + fz^2 = g^2 with least mean square.
% Geoffrey Vincent & Philipp Clausen
% 2015_07_08

% observations, all your observations  fx fy fz goes here.
disp('==============================');
fprintf('Il y a %.0f observations \n', length(obs))

% Approched value, you need non zeros if the model is non linear.
kx = init(1); 
bx = init(2);

ky  = init(3);
by  = init(4);
bxy = init(5);

kz  = init(6);
bz  = init(7);
bxz = init(8);
byz = init(9);

% %==========================================================================
% % Modele stochastique: matrice de covariance des observations
% sigmafx = 0.11;
% sigmafy = 0.12;
% sigmafz = 0.13;
% 
% Qll = eye(length(obs),length(obs));
% 
% for i = 1:lenght(obs)
%     
%     Qll(3*i-2, 3*i-2) = sigmafx^2;
%     Qll(3*i-1, 3*i-1) = sigmafy^2;
%     Qll(3*i,     3*i) = sigmafz^2;
%     
% end
% 
% P = inv(Qll);
P = eye(length(obs));
Qll = inv(P);
% %=========================================================================

% Start compensation iterative
% Prealocation
it = 0;
dx = ones(6,1);
A  = zeros(length(obs), 6);
L = zeros(length(obs), 1);
tic();
while max(abs(dx))> 0.00000001 && it < 50
   
   it = it + 1;

   for i = 1:length(obs)
        
        % Calcul des derivees partielle 
        [a1, a2, a3, a4, a5, a6, e1] = ...
            ModelFontionnel(kx, bx, ky, by, bxy, kz, bz, bxz, byz, obs(i,:));
        
        % Remplissage de la matrice Jacobienne avec les derivees partielles
        A(i, 1:6) = [a1, a2, a3, a4, a5, a6];

        % Remplissage de la matrice des ecarts apparents
        L(i, 1) = e1; % mm
   end

   % Covariance des parametres compenses Qxx
   Qxx = (A' * P * A)^-1;

   % Ajustement vector dx
   dx = Qxx * A' * P * L;
   
   bx  = bx  + dx(1);
   by  = by  + dx(2);
   bxy = bxy + dx(3);
   bz  = bz  + dx(4);
   bxz = bxz + dx(5);
   byz = byz + dx(6);
   
end

% Correlations des parametres compenses
for i = 1:size(A, 2)
    for j = 1:size(A, 2)
        corr_x(i, j) = Qxx(i,j)/sqrt(Qxx(i,i)*Qxx(j,j));
    end
end

% Create a colored plot of the correlation matrix
figure(99)
clf
imagesc(abs(corr_x));            
colormap((gray)); 
colorbar
drawnow

V = A*dx - L; % V, matrice des residus compenses
mqo = sqrt( V'*P*V / (size(A, 1)-size(A, 2)) );

Qllcomp = A * Qxx * A';
Qvv     = Qll - Qllcomp;
sigmax  = sqrt(diag(Qxx)) * mqo;
sigmaV  = sqrt(diag(Qvv)) * mqo;
quoloc  = V ./ sigmaV;

    fprintf('Results after %1.0f iterations in %0.2f s: \n', it, toc());
    fprintf('bx  = %.7f +/- %.5f \n', bx,  sigmax(1))
    fprintf('by  = %.7f +/- %.5f \n', by,  sigmax(2))
    fprintf('bxy = %.7f +/- %.5f \n', bxy, sigmax(3))
    fprintf('bz  = %.7f +/- %.5f \n', bz,  sigmax(4))
    fprintf('bxz = %.7f +/- %.5f \n', bxz, sigmax(5))
    fprintf('byz = %.7f +/- %.5f \n', byz, sigmax(6))
    disp('  ');
    
    param = [kx bx ky by bxy kz bz bxz byz];
end

function [a1, a2, a3, a4, a5, a6, e1] = ModelFontionnel(kx, bx, ky, by, bxy, kz, bz, bxz, byz, obs)
% Derivation des fonction utilisees dans la matrice jacobienne A (appelee aussi modele fonctionnel)
% La methode utilisee ici est la derivation numerique.

fx   = obs(1);
fy   = obs(2);
fz   = obs(3);
h    = 0.000000001; % Infinitesimal value
graw = computeg(kx, bx, ky, by, bxy, kz, bz, bxz, byz, fx, fy, fz);

a1 = (computeg(kx, bx+h, ky, by, bxy, kz, bz, bxz, byz, fx, fy, fz) - graw) / h;
a2 = (computeg(kx, bx, ky, by+h, bxy, kz, bz, bxz, byz, fx, fy, fz) - graw) / h;
a3 = (computeg(kx, bx, ky, by, bxy+h, kz, bz, bxz, byz, fx, fy, fz) - graw) / h;
a4 = (computeg(kx, bx, ky, by, bxy, kz, bz+h, bxz, byz, fx, fy, fz) - graw) / h;
a5 = (computeg(kx, bx, ky, by, bxy, kz, bz, bxz+h, byz, fx, fy, fz) - graw) / h;
a6 = (computeg(kx, bx, ky, by, bxy, kz, bz, bxz, byz+h, fx, fy, fz) - graw) / h;

% Calculs de l'ecart apparent e1
% gv = 9.80665;
gv = 9.8055; % Jan Skalouds gravity according to Lab 6 of Sensor Orientation
e1 = gv - graw;

end

function gcalc = computeg(kx, bx, ky, by, bxy, kz, bz, bxz, byz, fx, fy, fz)
% Compute g value with g = (fx^2 + fy^2 + fz^2)^0.5

    gcalc = sqrt( kx^2*(fx+bx)^2 + (ky*(fy+by) + bxy*kx*(fx+bx))^2 + ...
                 (kz*(fz+bz) + bxz*kx*(fx+bx) + byz*(ky*(fy+by) + bxy*kx*(fx+bx)))^2 );

end
