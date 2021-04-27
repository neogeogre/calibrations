%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = calib()
% Calibration of IMU accelerometers
% Denis Rouzaud, December 2007
% Projet SIE, Topo, EPFL
%---------------------------------
% Mabillard Romain
% 2013.03 : modified to allow the use of readimu.m instead of readFSAS.m 
%---------------------------------
% Philipp Clausen 
% 2015.02 : major modifications for gravity value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PREPARE AND CLEAR EVERYTHING
clear all
close all
fclose('all')
clc


%% HOW MUCH THE DATA SHALL BE DECORALATED?
%use 1 meas for [dms] meas
dms = 10;


%% FULL (=1) OR PARTIAL (=0) ESTIMATION 
% FULL: MISALIGNEMENT, SCALE FACTOR, BIAS
% PARTIAL: MISALIGNEMENT, BIAS
estimation_mode = 1;


%% DEFINE THE TYPE OF DATA
imuType=input('IMU type (IMAR,LN200, LN200IG,IXSEA, XSENS, NAVCHIP_INT, NAVCHIP_FLT) : ','s');


%% THE ACTUAL PROCESS

%tolerance [m/s2]
%used in checkdata
tolerance = .08;

%reference signals
% g = 9.7955229; %m/s2
g = 9.8055; %m/s2 according to Jan Skalouds Sensor Orientation Lab 6

%step max to fix time [s] in IMU files
%used in readdata
step_max = 0.1;


while(1)
    %Read index file
    [id_pos,time_start,time_duration,id_file,filename,x,y,z,stop,msg,fput,fputfile] = index();
    if stop == 1
        msgbox(msg,'Error','error');
        break
	end
	
    %Read data file(s)
    [fdata,fidpos,fidpost,stop,msg,freq] = readdata(imuType,id_pos,time_start,time_duration,id_file,filename,step_max);
    
    if stop == 1
        msgbox(msg,'Error','error');
        break
    end
    
    step=1/freq;
    
	checkon = menu('Would you like to check data ?','Yes','No');
	if checkon == 1
		[check,msg] = checkdata(fdata,tolerance,step);
		if check == 0
			msgbox(msg,'Error','error');
			pause(.6)
			cont = menu('Would you like to continue with these measures ?','Yes','No');
			if cont == 2
				break
			end
		end
	end

	ii = length(fdata);
	sigma = zeros(1,3);
	for ax = 1:3
		sigma(ax) = 0;
		for i=1:ii
			sigma(ax) = sigma(ax) + std(fdata{i}(:,ax));
		end 
		sigma(ax) = sigma(ax)/(ii*step);
	end
	sigma
	
	% assembling Measures
	lgx = [];
	lgy = [];
	lgz = [];
	for i = 1:length(fdata)
		lgx = [lgx; fdata{i}(:,1)];
		lgy = [lgy; fdata{i}(:,2)];
		lgz = [lgz; fdata{i}(:,3)];
	end
	l = [lgx lgy lgz];

	%selecting measures
	if dms ~= 1
		nms = length(l(:,1));
		l_less = zeros(floor(nms/dms),3);
		for j = 1:dms:nms
			l_less(1+(j-1)/dms,:) = l(j,:);
		end
		l = l_less;
		clear l_less
	end

	%transfrom m/s to m/s2 (integration)
%     l = l * freq;
    
    vec=[0
         0
         0
         0
         0
         0];

    %compensation
% % % % %     [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation(l,g,sigma,step,vec);
% % % % %  
% % % % % %     [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation_biasscale(l,g,sigma,step,vec);
% % % % % 
% % % % % %     it = length(Qxx);
% % % % % %     TL = [{'bx'},{'by'},{'bz'},{'Sx'},{'Sy'},{'Sz'},{'tetyz'},{'tetzx'},{'tetzy'}];
% % % % % % 
% % % % % %     imagesc(abs(Corr{it}))
% % % % % %     colormap gray
% % % % % %     axis square
% % % % % %     set(gca,'Xtick',1:9,'Ytick',1:9,'XtickLabel',TL,'YtickLabel',TL)
% % % % % %     colorbar
% % % % % %     title('Correlation matrix')
    if estimation_mode == 1
        [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation(l,g,sigma,step,vec);
    elseif estimation_mode == 0
        [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation_biasscale(l,g,sigma,step,vec);
    else
        error
    end
    
	output(xcomp,dx,Qxx,Corr,v,sigmapos,sigma,fput,fputfile);
	
    break
end