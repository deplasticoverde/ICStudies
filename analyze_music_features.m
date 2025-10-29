function [rhythmicity, melodicity, harmonicity] = analyze_music_features(audiofile)
% ============================================================
% ANALYZE_MUSIC_FEATURES
% Estima medidas simples de ritmicidad, melodicidad y armonicidad.
% 28-10-2025
% Uso:
%   [rhythmicity, melodicity, harmonicity] = analyze_music_features('archivo.wav')
%
% Descripción de índices:
%   - Ritmicidad: claridad del pulso en la envolvente temporal.
%   - Melodicidad: estabilidad espectral promedio.
%   - Armónicidad: regularidad del espaciado entre picos del espectro.
%
% Autor: Carolina Espinoza Oñate (LAÚD, Universidad de Chile)
% ============================================================

%% === 0. CARGAR Y PREPARAR AUDIO ===
[x, Fs] = audioread(audiofile);
if size(x, 2) > 1
    x = x(:,1); % convertir a mono
end

% Normalización RMS
x = x / sqrt(mean(x.^2));
L = length(x);
t = (0:L-1) / Fs;

fprintf('\nAnalizando: %s (Fs = %d Hz, duración = %.2f s)\n', audiofile, Fs, L/Fs);

%% 1. RITMICIDAD
% Basada en la autocorrelación de la envolvente de la señal.
env = abs(hilbert(x));       % envolvente analítica
env = env - mean(env);        % eliminar componente DC

[acor, lag] = xcorr(env, 'coeff');
acor = acor(lag >= 0);        % mitad positiva
lag  = lag(lag >= 0) / Fs;    % convertir a segundos

% Buscar primer pico significativo (excluyendo 0)
[pks, ~] = findpeaks(acor, ...
    'MinPeakProminence', 0.05, ...
    'MinPeakDistance', round(0.1 * Fs));

if isempty(pks)
    rhythmicity = 0;
else
    rhythmicity = max(pks);
end

fprintf('Ritmicidad (claridad de pulso): %.2f\n', rhythmicity);

%% 2. MELODICIDAD
% Basada en la estabilidad del centroide espectral (menor var = más tonal).
win  = hamming(round(0.05 * Fs)); % ventana de 50 ms
hop  = round(0.025 * Fs);         % salto de 25 ms
nfft = 2048;

nFrames   = floor((L - length(win)) / hop);
centroids = zeros(1, nFrames);

for i = 1:nFrames
    idx   = (1:length(win)) + (i - 1) * hop;
    frame = x(idx) .* win;
    X     = abs(fft(frame, nfft));
    X     = X(1:nfft/2);
    f     = (0:nfft/2 - 1) * (Fs / nfft);
    centroids(i) = sum(f(:) .* X(:)) / sum(X);
end

% Variación relativa del centroide
melodicity = 1 - std(centroids) / mean(centroids);
melodicity = max(0, min(1, melodicity)); % limitar entre 0–1

fprintf('Melodicidad (estabilidad espectral): %.2f\n', melodicity);

%% 3. ARMÓNICIDAD
% Basada en la regularidad del espaciado entre picos del espectro.
win = hamming(round(0.1 * Fs));   % ventana de 100 ms
frame = x(1:length(win)) .* win;  % primer frame representativo

X = abs(fft(frame, nfft));
X = X(1:nfft/2);
f = (0:nfft/2 - 1) * (Fs / nfft);
X = X / max(X);                   % normalización

% Detectar picos espectrales
[~, locs] = findpeaks(X, ...
    'MinPeakHeight', 0.05, ...
    'MinPeakDistance', 20);

if numel(locs) < 2
    harmonicity = 0;
else
    df = diff(f(locs));                % diferencias entre picos
    spread = std(df) / mean(df);       % dispersión relativa
    harmonicity = max(0, 1 - spread);  % 1 = picos regulares
end

fprintf('Armónicidad (estructura armónica): %.2f\n', harmonicity);

%% 4. GRÁFICOS
figure('Color', 'w', 'Name', audiofile);

% Señal temporal
subplot(3, 1, 1);
plot(t, x);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal temporal');

% Autocorrelación
subplot(3, 1, 2);
plot(lag, acor);
xlabel('Retardo (s)');
ylabel('Autocorrelación');
title('Autocorrelación de la envolvente');

% Barras resumen
subplot(3, 1, 3);
bar([rhythmicity, melodicity, harmonicity], 'FaceColor', [0.3 0.4 0.8]);
set(gca, 'XTickLabel', {'Ritmo', 'Melodía', 'Armonía'}, 'YLim', [0 1]);
ylabel('Índice (0–1)');
title('Medidas de musicalidad');

%% === 5. RESUMEN EN CONSOLA ===
fprintf('\nResumen:\n');
fprintf('   Ritmo:     %.2f\n', rhythmicity);
fprintf('   Melodía:   %.2f\n', melodicity);
fprintf('   Armonía:   %.2f\n\n', harmonicity);

end
