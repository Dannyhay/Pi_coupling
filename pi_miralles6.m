function [resa] = pi_miralles6(wints,winls,winss)
% wints=sumts;
% winls=sumls;
% winss=sumss;

% pi - Miralles et al.,2012
% compute constants
% data used here are air temperature, netradiation and evaporation

% get data size
[m,n,t] = size(wints);

% estimation of potential evapotranspiration
% A relatively mechanistic equation is was introduced by Preistly and Taylor (1972) as a simplified
% version of the more robust Penman (1948) and Penman-Monteith (1965) equations:
alpha = 1.26; % 1.26
rho_w = 1000; % kgm-3 density of water
cp = -1.013E3; % J kg-1 K-1
% gamma = 4.95E10-4; % psychrometric constant (kgm-3/degree cel)
% lambda_v = 2500; %kJ kg-1
% %  slope of the saturation vapor density curve(delta-(kg/degree cel))
% delta = (4098*(0.6108*exp(((17.27).*T)./(T+237.3))))./((T+237.3).^2);
% this is suggested in other literature


ff = m*n;
wints = reshape(wints,[],t);
winls = reshape(winls,[],t);
winss = reshape(winss,[],t);

% parpool(28)
parfor i = 1:ff
    
    T=squeeze(wints(i,:));Rn=squeeze(winls(i,:)); Ev =squeeze(winss(i,:));
    
    q_s    = 0.6108.*exp(17.3*T./(T+237.3));
	lambda = 1.91846*((T+273.15)/((T+273.15)-33.91)).^2;
    gamma  = cp./(0.622*lambda); % Brunt [1952]
    
    
    % delta=2508.3./(T+237.3).^2.*exp(17.3*T./(T+237.3));
	delta  = (q_s*4098)./((T+237.3).^2);
    
    % estimate net radiation from sensible and latent heat fluxes
    % Rn = SH + LH;
    RN(i,:)=Rn;
    
    % Potential evaporation
	PET = (alpha./(lambda*rho_w)).*((delta./(delta+gamma)).*Rn);
    pET(i,:) = PET;
    
    % computing pi
    %Pi = (1/std(T))*((cov(Ev,T)/std(Ev))-(cov(PET,T)/std(PET)));
    %res(i) = Pi(1,2);
	
	Pi = corr(T,Ev)-corr(T,PET);
	res(i) = Pi;
    
    disp([i]) 
    
end

resa = reshape(res,[m n]);
%delete(gcp('nocreate'))

 %mm = round(m/2);
 %Pi2 = cat(1,resa(mm:end,:),resa(1:mm-1,:));
 Pi2 = resa;
