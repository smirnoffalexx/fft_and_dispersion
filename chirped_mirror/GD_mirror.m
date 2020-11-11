%function []=GD_mirror()
tau = 5e-15; % pulse duration, s (FWHM)
dt = 2e-16;
t_max = 30e-14;
N = 2 * t_max / dt;
t = linspace (-t_max, t_max - dt, N);
% NFFT = 2^nextpow2 (length(t));
% t=0:dt:2*t_max-dt;
c = 3e8; % speed of light, m/s
E0 = 1; % field amplitude

dw = 2 * pi / ((length (t) - 1) * dt);
w_max = pi / dt;
w = linspace(-w_max, w_max - dw, N);
%w = 0:dw:2*w_max;


wavelength = 800e-9; % wavelenght, m
w0 = 2 * pi * c / wavelength; % angular frequency, central frequency rad/s
%lambda=675:15:875;
%GD=[1.5,0.25,-0.25,0,0.3,-0.25,0.35,-0.25,0.45,-0.3,0,48,-0.35,0.45,-0.45,0.45,-0.45]; 
%lambda=725:200/400:875-0.5; % new from 725 nm old from 675
%w=2*pi*c./lambda*10^9;

%%% Gaussian experimental Block
E_t = E0 * exp (-1/2 * (t/tau).^2) .* exp(1i * w0 * t);

subplot (3, 2, 1);
plot (t, abs (E_t).^2);
xlim ([-30e-15 30e-15]);
xlabel ('t, s');
ylabel ('Intensity, arb. units');
title ('Input pulse');
grid on;

% E_w = fft (E_t);
% plot (w, fftshift (abs (E_w).^2 / max (abs (E_w).^2)));
% xlim ([0.1e16 0.4e16]);
% grid on;
% 
% E_t1 = ifft (E_w);
% plot (t, abs (E_t1).^2);
% grid on;

N1 = 49;
lambda = linspace (725, 900, N1) * 1e-9;
GD_osc = 0.45 * sin (lambda / 4e-9);
% plot (lambda, GD_osc);

y0 = 144.16831;
A1 = - 2918.45365;
t1 = 2.15559e-7;
GD_exp = (y0 + A1 * exp (- lambda / t1));
% plot (lambda, GD_exp);

GD = GD_osc + GD_exp;
subplot (3, 2, 2);
plot (lambda, GD);
xlabel ('lambda, m');
ylabel ('GD, fs');
title ('GD approximation');
grid on;

GD_e = [17, 30, 43, 55, 65, 72, 80, 88, 93, 100]; % from graph
lambda_i = [675, 700, 725, 750, 775, 800, 825, 850, 875, 900]; % lambda from graph
subplot (3, 2, 3);
plot (lambda_i, GD_e);
xlabel ('lambda, nm');
ylabel ('GD, fs');
title ('GD from website');
grid on;

Phase_lambda = cumtrapz (lambda, 2e-15 * pi * c * GD .* (lambda).^(-2));
%Phase_lambda = - 0.45 * 4e-9 * cos (lambda / 4e-9) - A1 * t1 * exp (- lambda / t1);
subplot (3, 2, 4);
plot (lambda, Phase_lambda);
xlabel ('lambda, m');
ylabel ('Phase, rad');
title ('Phase from numerical calculation');
grid on;

%w = 2 * pi * c ./ lambda;

%plot (w, Phase_lambda);
%grid on;

Phase_w_detrend = detrend (Phase_lambda);
subplot (3, 2, 5);
plot (lambda, Phase_w_detrend);
xlabel ('lambda, m');
ylabel ('Phase, rad');
title ('Detrended Phase');
grid on;

w_before = zeros (1, 1700);
Phase_flip = fliplr (Phase_w_detrend);
w_after = zeros (1, length (w) - 1749);
Phase_w = [w_before , Phase_flip, w_after];
% plot (w, Phase_w);
% grid on;


% Phase_w = (-A * cos(pi/a * (w - w0 - xc)) * a/pi - A1 * x1 * exp (-(w - w0)/x1)) * 10^(-15); % Recalculated Phase taking into account w-w0
% Phase_detrend_w = detrend (Phase_w);
% E_w = fft(E_t);
% Ew = tau * sqrt (2 * pi) * exp (-((w - w0) * tau).^2 / 2) .* exp (1i * Phase_detrend_w); % Analytical spectrum calculation 
pulse = ifft (E_w .* exp(1i * Phase_w));
subplot (3, 2, 6);
plot (t, abs (pulse).^2 / max (abs (pulse).^2));
xlim ([-30e-15 30e-15]);
xlabel ('t, s');
ylabel ('Intesity, arb. units');
title ('Pulse after & before transformation');
grid on;
hold on;
plot (t, abs (E_t).^2);

%GD_e=[17,30,43,55,65,72,80,88,93, 100]; % from graph
%lamda_i=[675,700,725,750,775,800,825,850,875,900]; % lambda from graph 
