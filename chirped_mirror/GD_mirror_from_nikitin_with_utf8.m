clear
close all

tau = 15e-15; % pulse duration, s (FWHM)
 
t_max = 300e-15;
N = 256; % 256;
t = linspace (- t_max, t_max, N);
dt = t (2) - t (1);
c = 3e8; % speed of light, m/s
lambda0 = 800e-9; % wavelength, m

w = linspace (- pi / dt, pi / dt, N); % �������� ���� � ������. ������� ������ ����� ������ ��� fftshift
dw = w (2) - w (1);
w0 = 2 * pi * c / lambda0;               % angular frequency, central frequency rad/s

%%% Gaussian experimental Block
E0 = 1; % field amplitude
E_t = E0 * exp (- t.^2 / (2 * tau^2));
% E_w = fftshift(fft(E_t));
E_w = fft (E_t);
E_w_shifted = fftshift (E_w);
I_lambda = (abs (E_w_shifted)).^2 .* (w + w0).^2 / (2 * pi * c);    % �������� � ����� ����� � ������ �����������  ����� �����
lambda = 2 * pi * c ./ (w + w0);

% ��� � ���� ����� ����� �������� - �������� ��������� �������� ������� ��
% �������, ������������ - ���� ������������� ����������� 

GD_osc = 0.45 * sin (lambda / 4e-9);
% plot (lambda, GD_osc);

y0 = 144.16831;
A1 = - 2918.45365;
t1 = 2.15559e-7;
GD_exp = (y0 + A1 * exp (- lambda / t1));
% plot (lambda, GD_exp);
GD = GD_osc + GD_exp; % detrend (GD_osc + GD_exp);
phase = cumtrapz (lambda, 2e-15 * pi * c .* GD .* (lambda).^(-2));

% phase = 20 *(1e-15*w)+50*(1e-15*w).^2; % radians ������������ ������ �������
% phase = detrend(phase);   % ����� �������� �� ������� ��� �������������
% GD = diff(phase)./dw;       % ��������������� ����� �������� � ��������
% GD = [GD(1) GD];            % ��������� ����� ������� ����� ���� ��� ����� ������� �������

E2_w = E_w .* exp (1i * 1 * fftshift (detrend (phase))); % detrend (phase));
E2_t = ifft (E2_w);


% for Ti:Sa
B1 = 1.43134930;
B2 = 0.65054713; 
B3 = 5.3414021;
C1 = 5.2799261e-3 * 10^(-12);
C2 = 1.42382647e-2 * 10^(-12);
C3 = 325.017834 * 10^(-12);
n_lambda = (1 + B1 * lambda.^2 / (lambda.^2 - C1) + B2 * lambda.^2 / (lambda.^2 - C2) + B3 * lambda.^2 / (lambda.^2 - C3))^(1/2);
k_w = fliplr (n_lambda) * w / c;
L_crystal = 5e-3;
GD_crystal = L_crystal * diff (k_w)./dw;
GD_crystal = [GD_crystal(1) GD_crystal];
GD_sum = detrend (GD_crystal - 1 * fliplr (GD));

figure

 subplot (2,2,1);
 plot (lambda * 1e9, GD * 1e15);
 grid on;
 xlim ([600 1000]);
 xlabel ('nm');
 ylabel ('GD, fs');
 
 subplot (2,2,2);
 plot (lambda * 1e9, (I_lambda));
 xlim ([600 1000]);
 grid on;
 title ('C�����');
 xlabel ('nm');
 
 subplot(2,2,3);
 plot (t * 1e15, (abs (E_t)).^2);
 hold on;
 plot (t * 1e15, (abs (E2_t)).^2, 'r');
 grid on;
 xlabel ('fs');
 legend ('��', '�����');
 
 subplot (2,2,4);
 plot (w, GD_sum);
 %xlim ([(w0 - 5e14) (w0 + 5e14)]);
 grid on;
 title ('GD\_crystal - GD');
 xlabel ('w');