%% === 1. Setup iniziale ===
clear; clc; close all;

% Parametri soggetto
altezza = 1.74;       % Altezza (m)
massa = 73;           % Massa (kg)
cartella = 'OneDrive_1_05-05-2025';

fileNames = { ...
    'Arthur trial 1.xlsx', ...
    'Arthur trial 2.xlsx', ...
    'Arthur trial 3.xlsx', ...
    'Arthur trial 4.xlsx', ...
    'Arthur trial 5.xlsx' ...
};

nTrials = numel(fileNames);
trial_names = erase(fileNames, '.xlsx');
colori = lines(nTrials);

% Inizializza variabili di output
area_YZ_m2 = zeros(1, nTrials);
volume_3D_m3 = zeros(1, nTrials);
v_media = zeros(1, nTrials);
energia_proxy = zeros(1, nTrials);
energia_Joule = zeros(1, nTrials);

%% === 2. Inizializza figure ===
figure(1); hold on;
title('Traiettoria CoM nel piano YZ [m]');
xlabel('Z [m]'); ylabel('Y [m]'); axis equal; grid on;

figure(2); hold on;
title('Volume 3D del CoM [m³]');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); axis equal; grid on;

%% === 3. Loop sui file ===
for i = 1:nTrials
    filePath = fullfile(cartella, fileNames{i});
    opts = detectImportOptions(filePath);
    opts.DataRange = 'A12';
    data = readmatrix(filePath, opts);

    % === Tempo ===
    tempo = data(:,2);                  % tempo in secondi
    dt_vect = diff(tempo);              % dt tra frame
    T_tot = tempo(end) - tempo(1);

    % === Marker bacino ===
    PSIS_R = data(:, 12:14);
    PSIS_L = data(:, 15:17);
    ASIS_R = data(:, 18:20);
    ASIS_L = data(:, 21:23);

    % === Centro del bacino ===
    PelvisCenter = (PSIS_R + PSIS_L + ASIS_R + ASIS_L) / 4;

    % === Centro di massa dinamico ===
    CoM = PelvisCenter;
    CoM(:,2) = PelvisCenter(:,2) * 0.55;  % solo Y modificato secondo altezza

    % === Calcolo velocità lungo Z frame per frame ===
    dZ = diff(CoM(:,3));                       % spostamento lungo Z
    velocita_ist = abs(dZ ./ dt_vect);         % velocità istantanea (Z)

    v_media(i) = mean(velocita_ist);                   % velocità media Z
    energia_proxy(i) = mean(velocita_ist.^2);          % energia/massa
    energia_Joule(i) = 0.5 * massa * energia_proxy(i); % energia (J)

    % === Volume 3D ===
    try
        [K3, volume] = convhull(CoM(:,1), CoM(:,2), CoM(:,3));
        volume_3D_m3(i) = volume;
    catch
        volume_3D_m3(i) = NaN;
    end

    % === Area piano YZ ===
    try
        [~, areaYZ] = convhull(CoM(:,3), CoM(:,2));
        area_YZ_m2(i) = areaYZ;
    catch
        area_YZ_m2(i) = NaN;
    end

    % === Plot YZ ===
    figure(1);
    plot(CoM(:,3), CoM(:,2), '-', ...
        'Color', colori(i,:), 'LineWidth', 1.5, ...
        'DisplayName', ['Trial ', num2str(i)]);

    % === Plot volume 3D ===
    figure(2);
    scatter3(CoM(:,1), CoM(:,2), CoM(:,3), 10, colori(i,:), 'filled');
    try
        trisurf(K3, CoM(:,1), CoM(:,2), CoM(:,3), ...
            'FaceColor', colori(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end

figure(1); legend show;
figure(2); legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

%% === 4. Grafici a barre ===
figure;
subplot(2,2,1);
bar(area_YZ_m2); title('Area YZ [m²]');
ylabel('Area [m²]'); xlabel('Trial');
xticklabels({'1','2','3','4','5'}); grid on;

subplot(2,2,2);
bar(volume_3D_m3); title('Volume 3D [m³]');
ylabel('Volume [m³]'); xlabel('Trial');
xticklabels({'1','2','3','4','5'}); grid on;

subplot(2,2,3);
bar(v_media); title('Velocità media lungo Z [m/s]');
ylabel('Velocità Z [m/s]'); xlabel('Trial');
xticklabels({'1','2','3','4','5'}); grid on;

subplot(2,2,4);
bar(energia_Joule); title('Energia stimata (Z) [J]');
ylabel('Energia (J)'); xlabel('Trial');
xticklabels({'1','2','3','4','5'}); grid on;

%% === 5. Tabella finale ===
T = table(trial_names(:), ...
          area_YZ_m2(:), ...
          volume_3D_m3(:), ...
          v_media(:), ...
          energia_proxy(:), ...
          energia_Joule(:), ...
          'VariableNames', {'Trial', 'Area_YZ_m2', 'Volume_3D_m3', ...
                            'Velocita_Z_media_m_s', 'Energia_proxy_Z_m2_s2', 'Energia_Z_Joule'});

disp('✅ Risultati finali considerando solo la direzione Z, calcolati frame per frame:');
disp(T);
