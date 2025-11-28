clc; clear; close all;

%% 1) Parámetros del LFSR
m          = 7;                  % longitud del registro
taps       = [7 6 3 1];          % taps del polinomio
estado_ini = [1 0 0 0 0 0 0];        % estado inicial (≠ todo ceros)
L          = 1000000;              % longitud de la secuencia de bits

% Generar bits b[n] con el LFSR (USANDO LA FUNCIÓN)
[b, estado_fin] = LFSR_seq(m, taps, estado_ini, L);

fprintf('Longitud = %d, media de b[n] = %.4f\n', L, mean(b));
disp('Estado final del LFSR (solo info):');
disp(estado_fin);

%% 2) TRANSMISOR: Modulación de amplitud 0 -> -A, 1 -> +A
A = 1;                     % amplitud (como dice el enunciado)

% Mapeo: 0 -> -1, 1 -> +1
s = A * (2*b - 1);         % s[n] en {-1, +1}

% Potencia de la señal s[n]
Ps = mean(s.^2);
fprintf('Potencia de la señal s[n]: Ps = %.4f\n', Ps);

%% 3) CANAL AWGN: x[n] = s[n] + w[n]

% SNR en dB: 0, 2, 4, ..., 14
SNRdB_vec = 0:2:14;
numSNR    = length(SNRdB_vec);

sigma2_vec = zeros(1, numSNR);
x_all      = zeros(numSNR, L);   % cada fila: x[n] para un SNR dado

for i = 1:numSNR
    SNRdB = SNRdB_vec(i);

    % Varianza del ruido: sigma^2 = Ps / 10^(SNRdB/10)
    sigma2 = Ps / (10^(SNRdB/10));
    sigma2_vec(i) = sigma2;

    % Ruido blanco Gaussiano w[n] ~ N(0, sigma^2)
    w = sqrt(sigma2) * randn(1, L);

    % Señal recibida
    x = s + w;

    x_all(i, :) = x;

    fprintf('SNR = %2d dB -> sigma^2 = %.5f\n', SNRdB, sigma2);
end

%% 4) DEMODULADOR y cálculo de probabilidad de error

Pe = zeros(1, numSNR);   % vector para guardar la probabilidad de error

for i = 1:numSNR
    % Señal recibida para este SNR
    x = x_all(i, :);
    
    % Regla de decisión:
    % b_hat[n] = 1 si x[n] > 0
    % b_hat[n] = 0 si x[n] <= 0
    b_hat = x > 0;    % da vector lógico con 0/1
    
    % Número de errores
    num_err = sum(b_hat ~= b);
    
    % Probabilidad de error para este SNR
    Pe(i) = num_err / L;
    
    % Mostrar en pantalla
    fprintf('SNR = %2d dB | sigma^2 = %.5f | errores = %5d | Pe = %.5e\n', ...
            SNRdB_vec(i), sigma2_vec(i), num_err, Pe(i));
end

% Tabla de resultados en el workspace (opcional)

resultados = [SNRdB_vec(:), sigma2_vec(:), Pe(:)];

resultados_tbl = table( SNRdB_vec(:), sigma2_vec(:), Pe(:), ...
    'VariableNames', {'SNR_dB', 'sigma2', 'Pe'});


% Si usas MATLAB, también puedes crear una tabla bonita:
% resultados_tbl = table(SNRdB_vec(:), sigma2_vec(:), Pe(:), ...
%     'VariableNames', {'SNR_dB', 'sigma2', 'Pe'});

% Gráfica de Pe (escala logarítmica) vs SNR (dB)

figure;
semilogy(SNRdB_vec, Pe, '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Probabilidad de error P_e');
title('Probabilidad de error vs SNR para canal AWGN');

