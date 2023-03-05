function [Pi2, pET] = pi_miralles_BergPET(wints,winls,winss,winrs,winSFW,winevp,ALT)


% pi - Miralles et al.,2012
% PET following Berg & McColl (2021)
% compute constants

% get data size
[m,n,t] = size(wints);

% estimation of potential evapotranspiration
% A relatively mechanistic equation is was introduced by Preistly and Taylor (1972) as a simplified
% version of the more robust Penman (1948) and Penman-Monteith (1965) equations:
%alpha = 1.26; % 1.26
%rho_w = 1000; % kgm-3 density of water
%cp = -1.013E3; % J kg-1 K-1
% gamma = 4.95E10-4; % psychrometric constant (kgm-3/degree cel)
% lambda_v = 2500; %kJ kg-1
% %  slope of the saturation vapor density curve(delta-(kg/degree cel))
% delta = (4098*(0.6108*exp(((17.27).*T)./(T+237.3))))./((T+237.3).^2);
% this is suggested in other literature

%------------------------------------------------------
% This script was written to verify the pi computations based on the Penman equations in Pi_Miralles5 script
% This verification is written following Berg & McColl 2021
%------------------------------------------------------
ff = m*n;
wints  = reshape(wints,[],t);
winls  = reshape(winls,[],t);
winss  = reshape(winss,[],t);
winrs  = reshape(winrs,[],t);
winSFW = reshape(winSFW,[],t);
winevp = reshape(winevp,[],t);
ALT    = reshape(ALT,[]);


% parpool(28)
parfor i = 1:ff
    
	tas  = squeeze(wints(i,:));
	alt  = ALT(i);
	huss = squeeze(winrs(i,:));
	hfls = squeeze(winls(i,:));
	hfss = squeeze(winss(i,:));
	SFW  = squeeze(winSFW(i,:));
	Ev   = squeeze(winevp(i,:));
	
	
	%### Latent heat of vaporatization
	Lv = 2500.8-2.36*(tas-273.16)+0.0016.*(tas-273.16).^2-0.00006.*(tas-273.16).^3;
	%#in J/g. Multiply 1000/1000000=1/1000 to get MJ/kg
	Lv = Lv/1000; %#MJ/kg
	%Lv2 = Lv * 1000000 %# J/kg
%	###################################3

	ps=101.3*((293-0.0065*alt)/293).^5.26; 
	vpsat = 0.611 * exp( 17.27 * (tas-273.16) / (tas-273.16 + 237.3) );  %# in kPa
	delta = vpsat*4098/(((tas-273.16)+237.3)^2);   %# in kPa/C
	gamma = 0.665*10^-3* (ps/1000);   %# in kPa/C ****CHECK LATER****
	%#gamma <- 1.013e-3/(Lv*0.622) * (ps/1000) ## accounting for varying Lv
	
	w = huss/(1-huss);
	es =  w/(w+0.622)*(ps/1000); %# in kPa, so the VPD is in kPa, too. 
	Lv = 2500.8-2.36*(tas-273.16)+0.0016*(tas-273.16)^2-0.00006*(tas-273.16)^3 ;
	%#in J/g. Multiply 1000/1000000=1/1000 to get MJ/kg
	Lv = Lv/1000;
	%print("PET")
	
	% cp = 1.013E3; %#J.g-1.K-1 = 1013J.kg-1.K-1 ; rho is in kg.m-3; vpsaat is in Pa=kg.m-1.s-2, chU in m.s-1;
	%#numerator=J.kg-1.K-1*kg.m-3*Pa*m.s-1=Pa.K-1.J.m-2.s-1
	
	%############## Open Water Penman - fron Shuttleworth 1993 handbook, equation 4.2.30
	pETr = 1/(Lv)*( delta/(delta + gamma)*(hfls+hfss)*0.0864 + gamma/(gamma+delta)*(6.43*(1+0.536*SFW)*(vpsat-es)));

	
    %############# Computing Pi ###################   
    % computing pi
    % Pi = (1/std(T))*((cov(Ev,T)/std(Ev))-(cov(PET,T)/std(PET)));
	Pi = corr(Ev,T)-corr(PET,T);
    res(i) = Pi;
    
    disp([i]) 
    
end

Pi2 = reshape(res,[m n]);

%delete(gcp('nocreate'))

