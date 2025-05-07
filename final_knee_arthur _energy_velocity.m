clear all;
clc;

% === Parametri soggetto ===
height = 1.74;     % m
mass = 73;         % kg
g = 9.81;          % m/s²
cartella = 'OneDrive_1_05-05-2025';  % Cartella dove sono i file

fileNames = {
    'Arthur trial 1.xlsx', ...
    'Arthur trial 2.xlsx', ...
    'Arthur trial 3.xlsx', ...
    'Arthur trial 4.xlsx', ...
    'Arthur trial 5.xlsx'
};

nTrials = numel(fileNames);

% === Inizializzazione della tabella per i risultati ===
T_results = table('Size', [nTrials, 13], ...
    'VariableTypes', repmat({'double'}, 1, 13), ...
    'VariableNames', {'Volume_m3', 'E_cinetica_J', 'E_potenziale_J', 'E_totale_J', ...
                      'Work_X_J', 'Work_Y_J', 'Work_Z_J', ...
                      'Perc_Work_X', 'Perc_Work_Y', 'Perc_Work_Z', ...
                      'Frame_picco_Z', 'Tempo_picco_Z_s', 'Vmax_Z_m_s'});
T_results.Trial = fileNames';

% Loop per ogni file
for i = 1:nTrials
    fname = fileNames{i};
    fprintf('\n========= %s =========\n', fname);

    T = readtable(fullfile(cartella, fname));
    
    % Estrazione delle posizioni dei marker del ginocchio
    lateral_epicondyle = [T{:,36}, T{:,37}, T{:,38}];  % Right lateral epicondyle of femur (X, Y, Z)
    fibula_head = [T{:,39}, T{:,40}, T{:,41}];          % Right head of the fibula (X, Y, Z)
    medial_epicondyle = [T{:,41}, T{:,42}, T{:,43}];    % Right medial epicondyle of femur (X, Y, Z)

    % Numero di frame
    num_frames = size(lateral_epicondyle, 1);

    % Calcolare il centroide (CoM) del ginocchio
    CoM_knee = zeros(num_frames, 3);
    for f = 1:num_frames
        CoM_knee(f,:) = (lateral_epicondyle(f,:) + fibula_head(f,:) + medial_epicondyle(f,:)) / 3;
    end

    % Trova il punto di massima flessione (minimo lungo Z)
    [~, idx_max_Z] = min(CoM_knee(:,3));  % Trova il punto con il minimo di Z (massima flessione)

    % Fase di discesa (fino al minimo di Z)
    CoM_knee_down = CoM_knee(1:idx_max_Z, :);

    % Calcolare la velocità lungo Z per il tratto di discesa
    velocity_z_down = diff(CoM_knee_down(:,3)) * (1 / T.Time(2));  % Velocità lungo Z per il tratto di discesa

    % Calcolare la velocità massima lungo Z nel tratto di discesa
    [v_max_down, idx_max_v_z_down] = max(abs(velocity_z_down));  % Velocità massima (assoluta) lungo Z
    t_max_v_z_down = T.Time(idx_max_v_z_down);  % Tempo del picco di velocità lungo Z
    CoM_max_v_z_down = CoM_knee_down(idx_max_v_z_down, :);  % Posizione del CoM al picco di velocità

    % Calcolare l'energia cinetica al punto di massima velocità lungo Z
    E_kin_max_v_z_down = 0.5 * mass * v_max_down^2;

    % === Calcolare le velocità lungo X, Y, Z ===
    velocity = zeros(num_frames-1, 3);
    for j = 2:num_frames
        velocity(j-1,:) = (CoM_knee(j,:) - CoM_knee(j-1,:)) * (1 / T.Time(2));  % Calcolare velocità tra frame
    end

    % === Calcolare l'energia cinetica per X, Y, Z ===
    kinetic_energy_x = 0.5 * mass * velocity(:,1).^2;
    kinetic_energy_y = 0.5 * mass * velocity(:,2).^2;
    kinetic_energy_z = 0.5 * mass * velocity(:,3).^2;

    % === Calcolare l'energia potenziale in Y ===
    potential_energy_changes_y = zeros(num_frames-1, 1);
    for j = 2:num_frames
        height_change_y = CoM_knee(j,2) - CoM_knee(j-1,2);
        potential_energy_changes_y(j-1) = mass * g * height_change_y;
    end

    % Solo per i cambiamenti positivi di energia potenziale
    positive_potential_energy_y = potential_energy_changes_y;
    positive_potential_energy_y(positive_potential_energy_y < 0) = 0;

    % === Calcolare il lavoro meccanico ===
    work_x = sum(abs(diff(kinetic_energy_x)));
    work_y_kinetic = sum(abs(diff(kinetic_energy_y)));
    work_y_potential = sum(positive_potential_energy_y);
    work_y = work_y_kinetic + work_y_potential;
    work_z = sum(abs(diff(kinetic_energy_z)));

    % Lavoro totale
    total_work = work_x + work_y + work_z;

    % Percentuale di energia per direzione
    percent_work_x = (work_x / total_work) * 100;
    percent_work_y = (work_y / total_work) * 100;
    percent_work_z = (work_z / total_work) * 100;

    % === Salvataggio dei risultati in tabella ===
    T_results.Volume_m3(i) = prod(max(CoM_knee_down) - min(CoM_knee_down));  % Volume
    T_results.E_cinetica_J(i) = sum(kinetic_energy_x) + sum(kinetic_energy_y) + sum(kinetic_energy_z);  % Energia cinetica totale
    T_results.E_potenziale_J(i) = sum(positive_potential_energy_y);  % Energia potenziale
    T_results.E_totale_J(i) = T_results.E_cinetica_J(i) + T_results.E_potenziale_J(i);  % Energia totale
    T_results.Work_X_J(i) = work_x;
    T_results.Work_Y_J(i) = work_y;
    T_results.Work_Z_J(i) = work_z;
    T_results.Perc_Work_X(i) = percent_work_x;
    T_results.Perc_Work_Y(i) = percent_work_y;
    T_results.Perc_Work_Z(i) = percent_work_z;
    T_results.Frame_picco_Z(i) = idx_max_v_z_down;
    T_results.Tempo_picco_Z_s(i) = t_max_v_z_down;
    T_results.Vmax_Z_m_s(i) = v_max_down;

    % === Output dei risultati ===
    fprintf('Lavoro totale: %.2f J | Energia totale: %.2f J\n', total_work, T_results.E_totale_J(i));
    fprintf('Lavoro X: %.2f J | Lavoro Y: %.2f J | Lavoro Z: %.2f J\n', work_x, work_y, work_z);
    fprintf('Percentuale energia X: %.1f%% | Y: %.1f%% | Z: %.1f%%\n', percent_work_x, percent_work_y, percent_work_z);
    fprintf('Velocità massima Z: %.2f m/s | Energia cinetica massima Z: %.2f J\n', v_max_down, E_kin_max_v_z_down);

    % === Visualizzazione animata 3D ===
    figure('Name', ['Animazione 3D - ', fname]);
    hold on;
    grid on;
    axis equal;

    % Cambia le coordinate per gli assi come richiesto
    % Assegniamo le coordinate in modo che l'asse X sia lungo Z e l'asse Z lungo X
    plot3(CoM_knee_down(:,3), CoM_knee_down(:,2), CoM_knee_down(:,1), 'k-', 'LineWidth', 2);
    scatter3(CoM_knee_down(:,3), CoM_knee_down(:,2), CoM_knee_down(:,1), 30, 'r', 'filled');

    % Mostra il punto di velocità massima lungo Z
    scatter3(CoM_max_v_z_down(3), CoM_max_v_z_down(2), CoM_max_v_z_down(1), 100, 'g', 'filled');

    xlabel('Z (m)');
    ylabel('Y (m)');
    zlabel('X (m)');
    title('Traiettoria CoM del Ginocchio (assi modificati) con velocità massima');
    
    % Creo gli handle per la leggenda con i colori corretti
    h1 = plot3(NaN, NaN, NaN, 'k-', 'LineWidth', 2);
    h2 = scatter3(NaN, NaN, NaN, 30, 'r', 'filled');
    h3 = scatter3(NaN, NaN, NaN, 100, 'g', 'filled');
    legend([h1, h2, h3], {'Traiettoria CoM', 'Punti CoM', 'Punto di Velocità Massima'}, 'Location', 'Best');

    % Animazione
    h_dot = plot3(CoM_knee_down(1,3), CoM_knee_down(1,2), CoM_knee_down(1,1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    h_path = plot3([], [], [], 'r-', 'LineWidth', 1.5);

    for f = 1:idx_max_Z
        % Aggiorna la posizione del punto del CoM
        set(h_dot, 'XData', CoM_knee_down(f,3), 'YData', CoM_knee_down(f,2), 'ZData', CoM_knee_down(f,1));

        % Aggiorna la traiettoria
        x_data = CoM_knee_down(1:f, 3);
        y_data = CoM_knee_down(1:f, 2);
        z_data = CoM_knee_down(1:f, 1);
        set(h_path, 'XData', x_data, 'YData', y_data, 'ZData', z_data);

        % Pausa per l'animazione
        pause(0.01);
    end
end

% === Riepilogo finale ===
disp('✅ RIEPILOGO COMPLETO:');
disp(T_results);
%%
filename = 'Risultati_Analisi_ Right_Knee_Arthur.xlsx';
writetable(T_results, filename, 'Sheet', 'Risultati');
fprintf('✓ Risultati salvati nel file Excel: %s\n', filename);