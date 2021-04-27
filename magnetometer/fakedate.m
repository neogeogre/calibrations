clc
clear all
% close all

s = rng;
nbrObs = 1000;

% number_of_staic_data = 1;

mag_fake_true(1,1:3) = [0.47 0 0];
N = 1;

for i = 2:1:nbrObs
    
    r = -pi + 2*pi.*rand(N,1);
    p = -pi + 2*pi.*rand(N,1);
    y = -pi + 2*pi.*rand(N,1);
    
    rotation = [cos(p)*cos(y),                          cos(p)*sin(y), -sin(p); ...
        sin(r)*sin(p)*cos(y) - cos(r)*sin(y),   sin(r)*sin(p)*sin(y) + cos(r)*cos(y),   sin(r)*cos(p); ...
        cos(r)*sin(p)*cos(y) + sin(r)*sin(y),   cos(r)*sin(p)*sin(y) - sin(r)*cos(y),   cos(r)*cos(p)];
    
    mag_fake_true(i,1:3) = (rotation*mag_fake_true(1,1:3)')';
end

%
% T = [   1.0096   -0.1115   -0.0403; ...
%         0.1281    1.0386    0.0818; ...
%         -0.0473   -0.0213    0.9752];
%
% K = [   20.3557    0           0;...
%         0          19.7130     0; ...
%         0          0           40.6120];
%
% B = [   -0.1945;...
%         -0.2051; ...
%         0.076];

T = [2    0    0; ...
     0.5  10    0; ...
     0.1  0.2  14];

B = [ 1;...
      -0; ...
      0.0];

mag_fake_distorted = nan(size(mag_fake_true));

% for i = 1:1:nbrObs
%     mag_fake_distorted(i,:) = (inv(T)*inv(K)*T*(mag_fake_true(i,:)'+B))';
% end

for i = 1:1:nbrObs
    
    mag_fake_distorted(i,:) = (inv(T)*(mag_fake_true(i,:)'+B))' + randn(1,3)*0.001;
    
end


% mag_fake_distorted_static = nan(number_of_staic_data*length(mag_fake_distorted),3);
%
% for i = 1:1:length(mag_fake_distorted)
%
%     distorted_vector = (inv(T)*inv(K)*T*(mag_fake_true(i,:)'+B))';
%     index_start = i*number_of_staic_data-(number_of_staic_data-1);
%
%     for j = 0:1:(number_of_staic_data-1)
%         mag_fake_distorted_static(index_start+j:index_start+j,1:3) = distorted_vector+0.0005*randn(1,3);
%     end
%
% end

%figure(1)
% clf
% subplot(4,1,1)
% plot(mag_fake_true(:,1:3))
% subplot(4,1,2)
% plot(sqrt(mag_fake_true(:,1).^2+mag_fake_true(:,2).^2+mag_fake_true(:,3).^2))
% subplot(4,1,3:4)
% plot3(mag_fake_true(:,1), mag_fake_true(:,2), mag_fake_true(:,3), '.g')
% set(gca, 'DataAspectRatio', [1 1 1]);
% axis equal

%figure(2)
% clf
% subplot(4,1,1)
% plot(mag_fake_distorted(:,1:3))
% subplot(4,1,2)
% plot(sqrt(mag_fake_distorted(:,1).^2+mag_fake_distorted(:,2).^2+mag_fake_distorted(:,3).^2))
% subplot(4,1,3:4)
% plot3(mag_fake_distorted(:,1), mag_fake_distorted(:,2), mag_fake_distorted(:,3), '.r')
% set(gca, 'DataAspectRatio', [1 1 1]);
% axis equal

rng(s);


% addpath('T:\common\Software\INS\readmag');
% mag = readmag('T:\common\DATA\navchip\2016_03_14\16.03.14_0.mag', 'NAVCHIP_FLT');


mag = readmag('T:\common\DATA\navchip\2015_12_08_multiposition\Afternoon2\ext_0.mag', 'NAVCHIP_FLT');
%obs = mag(1:2:6000,2:4);

obs = mag_fake_distorted;

figure(3)
subplot(4,1,1)
plot(obs(:,1:3))
subplot(4,1,2)
plot(sqrt(obs(:,1).^2+obs(:,2).^2+obs(:,3).^2))
subplot(4,1,3:4)
plot3(obs(:,1), obs(:,2), obs(:,3), '.g')
set(gca, 'DataAspectRatio', [1 1 1]);
axis equal
axis equal



init = [1 0 1 0 0 1 mean(obs(:,1)) mean(obs(:,2)) mean(obs(:,3))]';
% init = [1 0 1 0 0 1 0 0 0]';

param = compensMagnetometerSimple9(obs, init);

%%
calibrater_matrice = [param(1,1)    0   0; ...
                      param(1,2)    param(1,3)  0; ...
                      param(1,4)    param(1,5)  param(1,6)];
calibrater_bias = [ param(1,7);...
                    param(1,8);...
                    param(1,9)];                  

calibrated_obs = nan(size(obs));                
                
for i = 1:1:length(obs)
    
    calibrated_obs(i,:) = (calibrater_matrice*(obs(i,:)')-calibrater_bias)';
    
end

figure(4)
subplot(4,1,1)
plot(calibrated_obs(:,1:3))
subplot(4,1,2)
plot(sqrt(calibrated_obs(:,1).^2+calibrated_obs(:,2).^2+calibrated_obs(:,3).^2))
subplot(4,1,3:4)
plot3(calibrated_obs(:,1), calibrated_obs(:,2), calibrated_obs(:,3), '.g')
set(gca, 'DataAspectRatio', [1 1 1]);
axis equal
axis equal


