% Levenberg-Marquardt example/test for acc calibration
% Philipp Clausen
% 2015_07_08
% ---------------

%==========================================================================
% prepare everything
format long G
clear all
close all
clc

choiceModel = questdlg('Modele de calibration :', 'Choix du modele de calibration', ...
    '6 parametres','9 parametres','9 parametres');

switch choiceModel
    case '6 parametres'
        calibModel = '6param';
    case '9 parametres'
        calibModel = '9param';
end

%==========================================================================
% generate reference data
n_attitudes = 300;
% g           = 9.80665;         % [m/s^2]
g = 9.8055; % Jan Skalouds gravity according to Lab 6 of Senor Orientation
freq        = 250;             % [Hz]
t_end       = 1/250*1;         % [s]
time_vector = 0:1/freq:t_end-1/freq;
nlength     = length(time_vector);
noise_level = 0.01;
acc.leveled = [ 0*ones(size(time_vector))' 0*ones(size(time_vector))' g*ones(size(time_vector))'; ];

acc.data = [];
attitude = [];
for i = 1:1:n_attitudes
    
    attitude = [attitude; rand(1,3)*2*pi];
    
    R_to_outlevel_the_leveled = computeR(attitude(i,1), attitude(i,2), attitude(i,3));
    
    data = (R_to_outlevel_the_leveled*acc.leveled')' + noise_level*(rand(nlength,3)-0.5);
    
    acc.data = [acc.data ; data];
    
end

%==========================================================================
% True value for generating fake data
kx = 1+1e-4;
bx = 2.00867e-02;

ky = 1+4e-4;
by = 1.55483e-01;
bxy = 1e-4;

kz = 1-3.2e-3;
bz =-1.62864e-02;
bxz = -1.5e-4;
byz = 4.2e-4;

acc.measurement(:,1) = acc.data(:,1) / kx - bx;
acc.measurement(:,2) = ( acc.data(:,2) - bxy*acc.data(:,1) )/ky - by;
acc.measurement(:,3) = ( acc.data(:,3) - bxz*acc.data(:,1) - byz*acc.data(:,2) )/kz - bz;

fprintf('True value : \n');
fprintf('kx  = %.7f  \n', kx)
fprintf('bx  = %.7f  \n', bx)
fprintf('ky  = %.7f  \n', ky)
fprintf('by  = %.7f  \n', by)
fprintf('bxy = %.7f  \n', bxy)
fprintf('kz  = %.7f  \n', kz)
fprintf('bz  = %.7f  \n', bz)
fprintf('bxz = %.7f  \n', bxz)
fprintf('byz = %.7f  \n', byz)
disp(' ');

% figure(1)
% clf
% subplot(2,3,1)
% hold on;
% plot(acc.data(:,1), '-b.')
% plot(acc.measurement(:,1), '-gx')
% hold off;
% grid on;
% subplot(2,3,2)
% hold on;
% plot(acc.data(:,2), '-b.')
% plot(acc.measurement(:,2), '-gx')
% hold off;
% grid on;
% subplot(2,3,3)
% hold on;
% plot(acc.data(:,3), '-b.')
% plot(acc.measurement(:,3), '-gx')
% hold off;
% grid on;
% subplot(2,3,4)
% hold on;
% plot(sqrt(acc.data(:,1).^2 + acc.data(:,2).^2 + acc.data(:,3).^2), 'b')
% plot(sqrt(acc.measurement(:,1).^2 + acc.measurement(:,2).^2 + acc.measurement(:,3).^2), 'g')
% hold off;
% grid on;

% save('measurement.mat', '-struct', 'acc', 'measurement')

%==========================================================================
% Approximate values for kx, bx, ky, by, bxy, kz, bz, bxz, byz
Approximate(1) = 1;
Approximate(2) = 0;

Approximate(3) = 1;
Approximate(4) = 0;
Approximate(5) = 0;

Approximate(6) = 1;
Approximate(7) = 0;
Approximate(8) = 0;
Approximate(9) = 0;

for k = 10:10:n_attitudes
    % Achtung! k must not be under 10! 10 is minimum for redondancy (9 parameters)
    
    if strcmp(calibModel,'6param')
        conpensValue = compens6param(acc.measurement(1:k,:), Approximate);
    else
        conpensValue = compens9param(acc.measurement(1:k,:), Approximate);
    end
    %     pause(0.1) 
end
% return
% time_vector = time_vector';
%
%
% figure(1)
% clf
% subplot(2,2,1)
% hold on;
% plot(acc.ref(:,1), acc.ref(:,2), 'b')
% plot(acc.leveled(:,1), acc.leveled(:,2), 'r')
%
% hold off;
% subplot(2,2,2)
% hold on;
% plot(acc.ref(:,1), acc.ref(:,3), 'b')
% plot(acc.leveled(:,1), acc.leveled(:,3), 'r')
%
% hold off;
% subplot(2,2,3)
% hold on;
% plot(acc.ref(:,1), acc.ref(:,4), 'b')
% plot(acc.leveled(:,1), acc.leveled(:,4), 'r')
%
% hold off;
% subplot(2,2,4)
% hold on;
% plot(acc.ref(:,1), sqrt(acc.ref(:,2).^2 + acc.ref(:,3).^2 + acc.ref(:,4).^2), 'b')
% plot(acc.leveled(:,1), sqrt(acc.leveled(:,2).^2 + acc.leveled(:,3).^2 + acc.leveled(:,4).^2), 'r')
%
% hold off;
% legend('noisy and turned', 'leveled')
%
% return
% %   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
% %   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
% %   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
%
% clear all
% hold off
% %randn('seed',0);	% specify a particular random sequence for msmnt error
%
% %epsPlots = 0;  formatPlot(epsPlots);		% 1: make .eps plots, 0: don't
%
% % *** For this demonstration example, simulate some artificial measurements by
% % *** adding random errors to the curve-fit equation.
%
% global	example_number
%
% example_number = 1;			  % which example to run.
%
% consts = [ ];                         % optional vector of constants
%
% Npnt = 100;				  % number of data points
%
% t = [1:Npnt]';				  % independent variable
% % true value of parameters ...
% if example_number == 1, p_true  = [ 20   10   1  50 ]'; end
% if example_number == 2, p_true  = [ 20  -24  30 -40 ]'; end
% if example_number == 3, p_true  = [  6   20   1   5 ]'; end
%
% y_dat = lm_func(t,p_true,consts);
% y_dat = y_dat + 0.1*randn(Npnt,1);	  % add random noise
%
% % range of values for basic paramter search
% p1 = 0.1*p_true(1):0.2*p_true(1):2*p_true(1);
% p2 = 0.1*p_true(2):0.2*p_true(2):2*p_true(2);
% p3 = 0.1*p_true(3):0.2*p_true(3):2*p_true(3);
% p4 = 0.1*p_true(4):0.2*p_true(4):2*p_true(4);
%
% % parameter search
% for ip2 = 1:length(p2);
%    for ip4 = 1:length(p4);
% 	pt = [ p_true(1)  p2(ip2) p_true(3) p4(ip4) ];
% 	delta_y = ( y_dat - lm_func(t,pt,consts) );
% 	X2(ip2,ip4) = (delta_y' * delta_y)/2;
%    end
% end
%
% figure(1); % ------------ plot shape of Chi-squared objective function
%  clf
%  mesh(p2,p4,log10(X2))
%   xlabel('p_2')
%   ylabel('p_4')
%   zlabel('log_{10}(\chi^2)')
% plotfile = ['lm_exampA',int2str(example_number),'.eps'];
% %if epsPlots, print(plotfile,'-solid','-color','-deps','-F:28'); end
%
% % *** Replace the lines above with a read-in of some
% % *** experimentally-measured data.
%
% % initial guess parameters  ...
% if example_number == 1, p_init  = [  5   2  0.2  10 ]';  end
% if example_number == 2, p_init  = [  4  -5  6    10 ]';  end
% if example_number == 3, p_init  = [ 10  50  5    5.6 ]';  end
%
% weight = Npnt/sqrt(y_dat'*y_dat);	  % sqrt of sum of data squared
%
% p_min = -10*abs(p_init);
% p_max =  10*abs(p_init);
%
% %{
% figure(3)
%  clf
%  plot(t,y_dat,'o');
%   xlabel('t')
%   ylabel('y(t)')
% %}
%
% % Algorithmic Parameters
% %         prnt MaxIter  eps1  eps2  epx3  eps4  lam0  lamUP lamDN UpdateType
%    opts = [  3,    100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2,    11,    9,        1 ];
%
% [p_fit,Chi_sq,sigma_p,sigma_y,corr,R2,cvg_hst] =  ...
% 	lm('lm_func',p_init,t,y_dat,weight,-0.01,p_min,p_max,consts,opts);
%
% y_fit = lm_func(t,p_fit,consts);
%
% disp('    initial    true       fit        sigma_p percent')
% disp(' -------------------------------------------------------')
% disp ([ p_init  p_true  p_fit sigma_p 100*abs(sigma_p./p_fit) ])
%
% n = length(p_fit);
%
% lm_plots ( n, t, y_dat, y_fit, sigma_y, cvg_hst );
%
