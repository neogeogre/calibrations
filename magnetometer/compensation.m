function [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation(l0,g,sigma,step,vec)




%Initial parameters
if sum(vec)~=0
    x=vec;
else
    x = zeros(9,1);
end
lcomp = l0;

%nms
nms = length(l0(:,1))



Qll = zeros(3,3);
for ax = 1:3
	Qll(ax,ax) = sigma(ax)^2;
end


P = inv(Qll);

barre = waitbar(0);

it = 1;
stop = 1;
while stop == 1
	[xcomp(it,:),dx(it,:),lcomp,Qxx{it},Corr{it},v(it,:,:),sigmapos(it)] = compensation_process(x,lcomp,nms,g,Qll,P,it,'',barre);
	if it >= 5
	    stop = menu('Re-calculate ?','Yes','No');
	end
	x = xcomp(it,:)';
	it = it + 1;
end
[xcomp(it,:),dx(it,:),lcomp,Qxx{it},Corr{it},v(it,:,:),sigmapos(it)] = compensation_process(xcomp(it-1,:)',l0,nms,g,Qll,P,it,'final ',barre);
close(barre);
end

function [xcomp,dx,lcomp,Qxx,Corr,v,sigmapos] = compensation_process(x,l,nms,g,Qll,P,it,msg,barre)

    %progress bar opening
	txt = ['Iteration ' int2str(it) '. Processing ' msg 'compensation.'];
    waitbar(0,barre,txt);

	
    N = zeros(9,9);
    u = zeros(9,1);

    for ms = 1:nms
        lgx = l(ms,1);
        lgy = l(ms,2);
        lgz = l(ms,3);

        [A,B,w] = AB(lgx,lgy,lgz,x,g);

        %Compensation
        m = 1/(B*Qll*B');
        N = N + A' * m * A;
        u = u + A'* m * w;
        
        if ms/10000 == round(ms/10000)
            waitbar(ms/(2*nms));
        end
    end

	
    dx = - inv(N)*u
	xcomp = x + dx
    Qxx = inv(N)

	Corr = zeros(9,9);
	for i = 1:9
		for j = 1:9
			Corr(i,j) = Qxx(i,j) / sqrt(Qxx(i,i)*Qxx(j,j));
		end
	end

	%measures residuals
	txt = ['Iteration ' int2str(it) '. Calculating ' msg 'measures residuals.'];
	waitbar(.5,barre,txt);
	v = zeros(nms,3);
	vtPv = 0;
	wk = 0;
	for ms=1:nms
        lgx = l(ms,1);
        lgy = l(ms,2);
        lgz = l(ms,3);
		[A,B,w] = AB(lgx,lgy,lgz,x,g);
		m = 1/(B*Qll*B');
		t = A * dx + w;
		v(ms,:) = - t * m * [B(1,1)*Qll(1,1) B(1,2)*Qll(2,2) B(1,3)*Qll(3,3)];        

		k = m * t;
		wk = wk + w * k;

		for o=1:3
			vtPv = vtPv + v(ms,o)*P(o,o)*v(ms,o);
		end

		if ms/10000 == round(ms/10000)
			waitbar(.5+ms/(2*nms));
		end
	end

	sigmapos = vtPv / (2*nms + 1)

	disp('-w''k == v''Pv');
	wk
	vtPv

	
	lcomp = l + v;
    %verif
    test = closegap(lcomp(:,1),lcomp(:,2),lcomp(:,3),xcomp,g);
    disp('Close gap |w|');
	norm(test)
	
end

function [w] = closegap(lgx,lgy,lgz,x,g)
        bgx = x(1);
        bgy = x(2);
        bgz = x(3);
        Sgx = x(4);
        Sgy = x(5);
        Sgz = x(6);
        tetgyz = x(7);
        tetgzx = x(8);
        tetgzy = x(9);

        
        %Pseudo realistic signal from measures
        gx = (lgx - bgx) / (1+Sgx);
        gy = (lgy - bgy + tetgyz * gx) / (1+Sgy);
        gz = (lgz - bgz - tetgzy * gx + tetgzx * gy) / (1+Sgz);

        %close gap
        w = gx.^2 + gy.^2 + gz.^2 - g^2;
end

function [A,B,w] = AB(lgx,lgy,lgz,x,g)
        bgx = x(1);
        bgy = x(2);
        bgz = x(3);
        Sgx = x(4);
        Sgy = x(5);
        Sgz = x(6);
        tetgyz = x(7);
        tetgzx = x(8);
        tetgzy = x(9);

        
        %Pseudo realistic signal from measures
        gx = (lgx - bgx) / (1+Sgx);
        gy = (lgy - bgy + tetgyz * gx) / (1+Sgy);
        gz = (lgz - bgz - tetgzy * gx + tetgzx * gy) / (1+Sgz);

        %close gap
        w = gx^2 + gy^2 + gz^2 - g^2;


        %initiates matrix size for compensation
        A = zeros(1,9);

        %Condition
        %f = gx2 + gy2 + gz2 - |g|2
        %df/dU = 2 *(gx.dgx/dU + gy.dgy/dU + gz.dz/dU)

        %df/dbgx
        A(1,1) = 2 * ( ...
                 gx * -1/(1+Sgx) + ...
                 gy/(1+Sgy) * -tetgyz/(1+Sgx) + ...
                 gz/(1+Sgz) * ( tetgzy/(1+Sgx) - tetgzx*tetgyz/((1+Sgy)*(1+Sgx)) ) ...
                 );
        %df/dbgy
        A(1,2) = 2 * ( ...
                 gy * -1/(1+Sgy) + ...
                 gz/(1+Sgz) * -tetgzx/(1+Sgy) ...
                 );
        %df/dbgz
        A(1,3) = 2 * ( ...
                 gz * -1/(1+Sgz) ...
                 );

        %df/dSgx
        A(1,4) = 2 * -(lgx-bgx)/((1+Sgx)^2) .* ( ...
                 gx + ...
                 gy * tetgyz/(1+Sgy) + ...
                 gz/(1+Sgz) * ( -tetgzy + tetgzx*tetgyz/(1+Sgy) ) ...
                 );

        %df/dSgy
        A(1,5) = 2 * -1/((1+Sgy)^2) * ( ...
                 gy * (lgy - bgy + tetgyz*(lgx-bgx)/(1+Sgx)) + ...
                 gz*tetgzx/(1+Sgz) .* ( lgy-bgy + tetgyz*(lgx-bgx)/(1+Sgx) ) ...
                 );

        %df/dSgz
        A(1,6) = 2 * -1/((1+Sgz)^2) * ( ...
                 gz * ( lgz-bgz - tetgzy*(lgx-bgx)/(1+Sgx) + tetgzx*((lgy-bgy)/(1+Sgy) + tetgyz/(1+Sgy) * (lgx-bgx)/(1+Sgx)) ) ...
                 );

        %df/dtetgyz
        A(1,7) = 2 * (lgx-bgx)/((1+Sgy)*(1+Sgx)) *( ...
                 gy + ...
                 gz * tetgzx/(1+Sgz) ...
                 );

        %df/dtetgzx
        A(1,8) = 2 * gz/(1+Sgz) * ((lgy-bgy)/(1+Sgy) + tetgyz/(1+Sgy) * (lgx-bgx)/(1+Sgx));

        %df/dtetgzy
        A(1,9) = 2 * gz/(1+Sgz) * -(lgx-bgx)/(1+Sgx);


        %Parameters
        %f = gx2 + gy2 + gz2 - |g|2
        %df/dlU = 2 *(gx.dgx/dlU + gy.dgy/dlU + gz.dz/dlU)


        B = zeros(1,3);
        %df/dlgx
        B(1,1) = 2 * ( ...
                 gx /(1+Sgx) + ...
                 gy/(1+Sgy) * tetgyz/(1+Sgx) + ...
                 gz/(1+Sgz) * ( -tetgzy/(1+Sgx) + tetgzx*tetgyz/((1+Sgy)*(1+Sgx)) ) ...
                 );
        %df/dlgy
        B(1,2) = 2 * ( ...
                 gy /(1+Sgy) + ...
                 gz/(1+Sgz) * tetgzx/(1+Sgy) ...
                 );
        %df/dlgz
        B(1,3) = 2 * ( ...
                 gz /(1+Sgz) ...
                 );

end






