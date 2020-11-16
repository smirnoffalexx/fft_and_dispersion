clear
close all

tau = 15e-15; % pulse duration, s (FWHM)
 
t_max = 200e-15;
N = 256;
t = linspace (- t_max, t_max, N);
dt = t (2)- t (1);
c = 3e8; % speed of light, m/s
lambda0 = 800e-9; % wavelength, m

w = linspace (- pi / dt, pi / dt, N); % Выбираем ноль в центре. Поэтому дальше везде помним про fftshift
dw = w (2) - w (1);
w0 = 2 * pi * c / lambda0;               % angular frequency, central frequency rad/s

%%% Gaussian experimental Block
E0 = 1; % field amplitude
E_t = E0 * exp (-t.^2/(2*tau^2));
% E_w = fftshift(fft(E_t));
E_w = fft (E_t);
E_w_shifted = fftshift (E_w);
I_lambda = (abs (E_w_shifted)).^2 .* (w + w0).^2 / (2 * pi * c);    % пересчет в длину волны с учетом центральной  длины волны
lambda = 2 * pi * c ./ (w + w0);

% Вот с этой фазой можно поиграть - линейное слагаемое сдвигает импульс по
% времени, квадратичное - дает дисперсионное расплывание 

GD_osc = 0.45 * sin (lambda / 4e-9);
% plot (lambda, GD_osc);

y0 = 144.16831;
A1 = - 2918.45365;
t1 = 2.15559e-7;
GD_exp = (y0 + A1 * exp (- lambda / t1));
% plot (lambda, GD_exp);
GD = GD_osc + GD_exp; % detrend (GD_osc + GD_exp);
phase = cumtrapz (lambda, 2e-15 * pi * c .* GD .* (lambda).^(-2));

% phase = 20 *(1e-15*w)+50*(1e-15*w).^2; % radians относительно центра спектра
% phase = detrend(phase);   % общая задержка по времени нам малоинтересна
% GD = diff(phase)./dw;       % соответствующая групп задержка в секундах
% GD = [GD(1) GD];            % Подгоняем длину вектора чтобы было как длина вектора частоты

E2_w = E_w .* exp (1i * fftshift (detrend (phase))); % detrend (phase));
E2_t = ifft (E2_w);

figure
 subplot (3,1,1);
 plot (lambda * 1e9, GD * 1e15);
 grid on;
 xlim ([600 1000]);
 xlabel ('nm');
 ylabel ('GD, fs');
 subplot (3,1,2);
 plot (lambda * 1e9, (I_lambda));
 xlim ([600 1000]);
 grid on;
 title ('Cпектр');
 xlabel ('nm');
 subplot(3,1,3);
 plot (t * 1e15, (abs (E_t)).^2);
 hold on;
 plot (t * 1e15, (abs (E2_t)).^2, 'r');
 grid on;
 xlabel ('fs');
 legend ('До', 'После');
 