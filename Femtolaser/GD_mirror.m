%function []=GD_mirror()
tau=10*10^(-15); % pulse duration, s (FWHM)
t=-30*10^(-15):2*10^(-16):(30-0.2)*10^(-15); % e
c=3*10^8; % speed of light, m/s
E0=1; % field amplitude

wavelength=800*10^(-9); % wavelenght, m
w0=2*pi*c/wavelength; % angular frequency, central frequency rad/s
%lambda=675:15:875;
%GD=[1.5,0.25,-0.25,0,0.3,-0.25,0.35,-0.25,0.45,-0.3,0,48,-0.35,0.45,-0.45,0.45,-0.45]; 
lambda=725:200/400:875-0.5; % new from 725 nm old from 675
w=2*pi*c./lambda*10^9;

%%% Gaussian experimental Block
E_t=E0*exp(-1/2*(t/tau).^2);%.*exp(1i*w0*t);
plot(t,abs(E_t).^2);
E_w=fftshift(E_t);
plot(lambda,abs(E_w).^2);
E_t1=ifft(E_w);
plot(t,E_t1);
%%%

GD = 0.45 * sin(2*pi*c/w/4);
GD_detrend = 0.45*sin(lambda/4); % new GD (without trend) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot(lambda,GD_detrend);
Phase_lambda = - 0.45*cos(lambda/4); % Phase (lambda)

variables = {'№ lambda w GD_detrend'};
points = 1:length(lambda);
data = [points' lambda' w' GD_detrend'];
sw = ['' sprintf('%s', variables{:}) sprintf('\n')];
sw = [sw sprintf('%2d %.2f %.2f %.2f\n', data')];
fid = fopen('GD.txt', 'wt'); % открыли файл для записи
fprintf(fid, '%s', sw);
fclose(fid);
type('GD.txt');


y0 = 4.28748*10^-5;
xc=-0.0062;
a=12.56648;
A = 0.4501;
GD_w = y0+A*sin(pi/a*(w-xc));
Phase_w = -A*cos(pi/a*(w-w0-xc))*a/pi; % recalculated Phase taking into account w-w0
E_w = fft(E_t);
Ew = tau*sqrt(2*pi)*exp(-((w-w0)*tau).^2/2).*exp(1i*Phase_w); 
pulse = ifft(E_w.*exp(1i*Phase_w));
max_pulse = max(pulse);
pulse_norm = pulse/max_pulse;
plot(t,abs(pulse_norm).^2);


GD_e=[17,30,43,55,65,72,80,88,93]; % from graph
lamda_i=[675,700,725,750,775,800,825,850,875]; % lamda from graph 

% y0=228.46;
% A=-29.47;
% x=-1.41*10^15;
% GD_t=(y0+A*exp(-(w-w0)/x)+A*cos((w-w0)/x)); % *10^15 так как фс
% phase=(-A*x*exp(-(w-w0)/x)+0*A/2*x*sin((w-w0)/x))*10^(-15); % *10^15 так как фс
% 
% fig1=plot(w,GD_e);
% %for ii=1:lenght(w)
%     Ew=tau*sqrt(2*pi)*exp(-((w-w0)*tau).^2/2).*exp(1i*phase);
% %end
% pulse= ifft(Ew);
% max=max(pulse);
% pulse_norm=pulse/max;
% fig2=plot(t,pulse_norm);
% %plot(w,exp(-w.^2*tau^2/2)*exp(phase))
% %end

