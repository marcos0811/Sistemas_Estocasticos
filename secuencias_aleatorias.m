Nb = 10000; % numero de bits
Tb = 1; % tiempo de bit
fs = 100; % frecuencia de muestreo
Ts = 1/fs; % pperiodo de muestreo
Ns = Tb*fs; % muestras por bit

% generacion de la se;al
bk = randi([0 1], 1, Nb); % bits equiprobables
ak = 2*bk - 1;
s = repelem(ak, Ns); % se;al PCM
t = (0:length(s)-1)*Ts;  

figure();
plot(t(1:10*Ns), s(1:10*Ns), 'LineWidth', 1.5);
xlabel('Tiempo (s)'); ylabel('Amplitud');
title('Segmento de la señal PCM Polar NRZ');
grid on; ylim([-1.5 1.5]);

% ergodicidad, la media y autocorrelacion
media_temporal = mean(s);
fprintf('Media Temporal Estimada: %.5f\n', media_temporal);

% convergencia de la media
media_acumulada = cumsum(s) ./ (1:length(s));
figure();
plot(media_acumulada, 'LineWidth', 1.2);
yline(0,'r--','Media Teórica');
xlabel('Muestras'); ylabel('Media');
title('Convergencia de la media temporal'); grid on;

% autocorrelacion
[Rss, lags] = xcorr(s, 'biased');
tau = lags * Ts;
R_teorica = (1 - abs(tau)/Tb);
R_teorica(abs(tau) > Tb) = 0;

figure();
plot(tau, Rss, 'b', 'LineWidth', 1.2); hold on;
plot(tau, R_teorica, 'r--', 'LineWidth', 1.5);
xlabel('\tau (s)'); ylabel('R_{SS}(\tau)');
title('Autocorrelación: Temporal vs Teórica');
legend('Simulada','Teórica'); grid on; xlim([-2*Tb 2*Tb]);

% PDF empírica
figure();
histogram(s, 'Normalization','pdf');
xlabel('Amplitud'); ylabel('Densidad de probabilidad');
title('Histograma de la señal PCM Polar NRZ');
grid on;

% varias realizaciones 
figure();
hold on;
for k = 1:4
    bk_aux = randi([0 1], 1, 50);
    ak_aux = 2*bk_aux - 1;
    s_aux = repelem(ak_aux, Ns);
    t_aux = (0:length(s_aux)-1)*Ts;
    plot(t_aux, s_aux);
end
xlabel('Tiempo (s)'); ylabel('Amplitud');
title('Diferentes realizaciones del proceso PCM Polar NRZ');
grid on; ylim([-1.5 1.5]);

% estimacion de PSD
L = length(s);
f = (-L/2 : L/2-1)*(fs/L);

% PSD Teórica
arg = f*Tb;
Sinc_m = sin(pi*arg + eps) ./ (pi*arg + eps); 
PSD_teorica = Tb * (Sinc_m).^2;

% metodo indirecto
N_fft = length(Rss);
f_ind = (-N_fft/2 : N_fft/2-1)*(fs/N_fft);
PSD_indirecta = real(fftshift(fft(Rss))) * Ts; 

% metodo directo
S_f = fftshift(fft(s));
PSD_directa = (abs(S_f).^2 * Ts) / length(s);

% comparacion de PSDs normalizadas
PSD_dir_n = PSD_directa / max(PSD_directa);
PSD_ind_n = PSD_indirecta / max(PSD_indirecta);
PSD_teo_n = PSD_teorica / max(PSD_teorica);

plot(f, 10*log10(PSD_dir_n), 'Color',[0.7 0.7 0.7]); hold on;
plot(f_ind, 10*log10(PSD_ind_n), 'b','LineWidth',1.2);
plot(f, 10*log10(PSD_teo_n), 'r--','LineWidth',1.5);
xlabel('Frecuencia (Hz)'); ylabel('PSD Normalizada (dB)');
title('Comparación de PSD (Ajustadas a 0 dB)');
legend('Directo','Indirecto','Teórica');
grid on; xlim([-4/Tb 4/Tb]); ylim([-60 5]);

% PSD Directa vs Promediada 
PSD_prom = movmean(PSD_directa, 200);
figure();
plot(f, 10*log10(PSD_directa), 'Color',[0.8 0.8 0.8]); hold on;
plot(f, 10*log10(PSD_prom), 'b','LineWidth',1.3);
plot(f, 10*log10(PSD_teorica),'r--','LineWidth',1.5);
xlabel('Frecuencia (Hz)'); ylabel('PSD (dB/Hz)');
title('Efecto del Promediado en la Estimación Espectral');
legend('Directa','Promediada','Teórica');
grid on; xlim([-4/Tb 4/Tb]); ylim([-50 5]);

% PSD terica lineal
figure();
plot(f, PSD_teorica, 'r','LineWidth',1.5);
xlabel('Frecuencia (Hz)'); ylabel('PSD');
title('PSD Teórica (Escala Lineal)');
grid on; xlim([-4/Tb 4/Tb]);

% ancho de banda
BW_Nulo = 1/Tb;

%BW Estimado 95 porciento de la potencia
idx_pos = f >= 0;
f_pos = f(idx_pos);
PSD_pos = PSD_directa(idx_pos); 

Energia_Acum = cumsum(PSD_pos);
Energia_Total = Energia_Acum(end);

idx_95 = find(Energia_Acum >= 0.95 * Energia_Total, 1);
BW_95 = f_pos(idx_95);

fprintf('--- Ancho de Banda ---\n');
fprintf('BW Teórico (Primer Nulo): %.2f Hz\n', BW_Nulo);
fprintf('BW Estimado (95%% Potencia): %.2f Hz\n', BW_95);
