clear
close all

tau = 15e-15; % pulse duration, s (FWHM)

t_max = 1 * 330e-15; % 330 * 2 && N * 2 % for t_max = 300e-15 && N = 451 - limit for w_max = pi /dt
N = 1 * 450; % 256;
N1 = 1 * 256;
t = linspace (- t_max, t_max, N); 
dt = t (2) - t (1); %%%%%%%%%%%%%%%%%%%%%%%%%%%% bag was there w_max = pi/dt, if w_max > w0 then - w_max + w0 < 0, so it is bad, dt must always be more than pi/w0
c = 3e8; % speed of light, m/s
lambda0 = 800e-9; % wavelength, m

t1 = linspace (- t_max, t_max, N1); 
dt1 = t1 (2) - t1 (1);
w1 = linspace (- pi / dt1, pi / dt1, N1); % ¬ыбираем ноль в центре. ѕоэтому дальше везде помним про fftshift
dw1 = w1 (2) - w1 (1);
w = linspace (- pi / dt, pi / dt, N); % ¬ыбираем ноль в центре. ѕоэтому дальше везде помним про fftshift
dw = w (2) - w (1);
w0 = 2 * pi * c / lambda0;               % angular frequency, central frequency rad/s

%%% Gaussian experimental Block
E0 = 1; % field amplitude
E_t = E0 * exp (- t.^2 / (2 * tau^2));
% E_w = fftshift(fft(E_t));
E_w = fft (E_t);
E_w_shifted = fftshift (E_w);
I_lambda = (abs (E_w_shifted)).^2 .* (w + w0).^2 / (2 * pi * c);    % пересчет в длину волны с учетом центральной  длины волны
lambda = 2 * pi * c ./ (w + w0); % + lambda0;% * ones (1, N); %(w + w0);
lambda1 = 2 * pi * c ./ (w1 + w0);
E_t1 = E0 * exp (- t1.^2 / (2 * tau^2));
% E_w = fftshift(fft(E_t));
E_w1 = fft (E_t1);
E_w_shifted1 = fftshift (E_w1);
I_lambda1 = (abs (E_w_shifted1)).^2 .* (w1 + w0).^2 / (2 * pi * c);


% ¬от с этой фазой можно поиграть - линейное слагаемое сдвигает импульс по
% времени, квадратичное - дает дисперсионное расплывание 
GD_osc1 = 0.45 * sin (lambda1 / 4e-9);
GD_osc = 0.45 * sin (lambda / 4e-9);
% plot (lambda, GD_osc);

y0 = 144.16831;
A1 = - 2918.45365;
ta = 2.15559e-7;
GD_exp = (y0 + A1 * exp (- lambda / ta));
GD_exp1 = (y0 + A1 * exp (- lambda1 / ta));
% plot (lambda, GD_exp);
GD = GD_osc + GD_exp; % detrend (GD_osc + GD_exp);
GD1 = GD_osc1 + GD_exp1;
phase = cumtrapz (lambda, GD .* (w + w0).^2 / (2e15 * pi * c));
phase1 = cumtrapz (lambda1, GD1 .* (w1 + w0).^2 / (2e15 * pi * c));

% phase = 20 *(1e-15*w)+50*(1e-15*w).^2; % radians относительно центра спектра
% phase = detrend(phase);   % обща€ задержка по времени нам малоинтересна
% GD = diff(phase)./dw;       % соответствующа€ групп задержка в секундах
% GD = [GD(1) GD];            % ѕодгон€ем длину вектора чтобы было как длина вектора частоты

E2_w = E_w .* exp (1i * 1 * fftshift (detrend (phase))); % detrend (phase));
E2_t = ifft (E2_w);


% for Ti:Sa
B1 = 1.43134930;
B2 = 0.65054713; 
B3 = 5.3414021;
C1 = 5.2799261e-3 * 10^(-12);
C2 = 1.42382647e-2 * 10^(-12);
C3 = 325.017834 * 10^(-12);
n_lambda = (1 + B1 * lambda.^2 ./ (lambda.^2 - C1) + B2 * lambda.^2 ./ (lambda.^2 - C2) + B3 * lambda.^2 ./ (lambda.^2 - C3)).^(1/2); % Sellmeier formula
% k_w = n_lambda * 2 * pi ./ lambda; % k_w = fliplr (n_lambda) .* w / c; % k_w = 2 * pi * n_lambda ./ lambda; %
L_crystal = 5e-3;
dlambda = lambda (2) - lambda (1);
v = c ./ n_lambda; % phase velocity
u = v (2 : N) - lambda (2 : N) .* diff (v) / dlambda; % group velocity
GD_crystal = L_crystal ./ u;
% GD_crystal = L_crystal * diff (k_w) / (lambda (2) - lambda (1)) .* (lambda (2 : N)).^2 / (2 * pi * c); % L_crystal * diff (k_w)./ dlambda; % GD_crystal = 
GD_crystal = [GD_crystal(1) GD_crystal];
GD_sum =  (GD_crystal * 1e15); % - 2 * GD); % GD.^6); % detrend (GD_crystal - 2 * fliplr (GD));

figure

    subplot (2,2,1);
    plot (lambda * 1e9, GD); % , lambda1 * 1e9, GD1); %* 1e15);
    grid on;
    xlim ([700 900]);
    xlabel ('lambda, nm');
    ylabel ('GD, fs');

    subplot (2,2,2);
    plot (lambda * 1e9, (I_lambda) / max (I_lambda), lambda1 * 1e9, (I_lambda1) / max (I_lambda1), 'r');
    % plot (lambda * 1e9, (I_lambda), lambda1 * 1e9, (I_lambda1), 'r');
    % plot (w, (abs (fftshift (E_w))).^2 / max ((abs (fftshift (E_w))).^2) * dt, w1, (abs (fftshift (E_w1))).^2 / max (abs (fftshift (E_w1)).^2) * dt, 'r');
    % plot (w, (abs (fftshift (E_w * 2 / N))).^2, w1, (abs (fftshift (E_w1 * 2 / N1))).^2, 'r');
    xlim ([700 900]);
    grid on;
    title ('Spectrum');
    xlabel ('lambda, nm');

    subplot(2,2,3);
    plot (t * 1e15, (abs (E_t)).^2);
    hold on;
    plot (t * 1e15, (abs (E2_t)).^2, 'r');
    grid on;
    xlabel ('t, fs');
    legend ('Before', 'After');

    subplot (2,2,4);
    % plot (w + w0, GD_sum); % plot (lambda * 1e9, GD_sum);
    plot (lambda * 1e9, GD_sum);
    xlim ([700 900]);
    % xlim ([(w0 - 0.26e15) (w0 + 0.34e15)]); % xlim ([700 900]); % xlim ([((lambda0 - 100e-9) * 1e9)  ((lambda0 + 100e-9) * 1e9)]);
    grid on;
    title ('GD\_crystal - GD');
    xlabel ('lambda, nm');
    
% figure % for checking phase and E_t point number independency
%     subplot (2, 1, 1);
%     plot (lambda * 1e9, phase, lambda1 * 1e9, phase1, 'r');
%     xlim ([700 900]);
%     xlabel ('lambda, nm');
%     ylabel ('phase');
%     
%     subplot (2, 1, 2);
%     plot (t * 1e15, E_t, t1 * 1e15, E_t1, 'r');
%     xlim ([-350 350]);