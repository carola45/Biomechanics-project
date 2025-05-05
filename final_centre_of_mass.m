%% centre of mass - ho fatto per un trial alla volta, ho visto che tu li gestivi tutti insieme io non sapevo farlo
clear all
clc

% Parametri soggetto
height = 1.74;       % Altezza (m)
mass = 73;           % Massa (kg)

filename = 'Arthur trial 2.xlsx';
T = readtable(filename);


% traiettoria del CoM
% volume che comprende la traiettoria del CoM
% baricentro come baricentro asis e pesis -- da definire dinamico frame per
% frame
% voglio calcolare energia, ma solo nel piano x,y perchè è quello spreco di
% energia che non serve 

% extract markers positions
ASIS_R = [T.RightASIS_X, T.RightASIS_Y, T.RightASIS_Z];
ASIS_L = [T.LeftASIS_X, T.LeftASIS_Y, T.LeftASIS_Z];
PSIS_R = [T.RightPESIS_X, T.RightPESIS_Y, T.RightPESIS_Z];
PSIS_L = [T.LeftPESIS_X, T.LeftPESIS_Y, T.LeftPESIS_Z];

% Number of frames
num_frames = size(ASIS_R, 1);

%% Calculate centre of mass as the pelvis centre
CoM = zeros(num_frames, 3);

for i = 1:num_frames
    % Calculate midpoint of all four markers for each frame
    CoM(i,:) = mean([ASIS_R(i,:); ASIS_L(i,:); PSIS_R(i,:); PSIS_L(i,:)], 1);
end

%% Calculate volume around CoM trajectory

% Find min and max coordinates to define bounding box
min_x = min(CoM(:,1));
max_x = max(CoM(:,1));
min_y = min(CoM(:,2));
max_y = max(CoM(:,2));
min_z = min(CoM(:,3));
max_z = max(CoM(:,3));

% Calculate dimensions
x_range = max_x - min_x;
y_range = max_y - min_y;
z_range = max_z - min_z;

% Calculate volume
com_volume = x_range * y_range * z_range;

% Calculate path length
com_path_length = 0;
for i = 2:num_frames
    segment_length = norm(CoM(i,:) - CoM(i-1,:));
    com_path_length = com_path_length + segment_length;
end


%% figures

% Create figure
figure('Position', [100, 100, 1200, 800]);

% Plot markers and CoM trajectory
subplot(3,2,[1,3,5]);
hold on;
plot3(ASIS_R(:,1), ASIS_R(:,2), ASIS_R(:,3), 'r.', 'MarkerSize', 5);
plot3(ASIS_L(:,1), ASIS_L(:,2), ASIS_L(:,3), 'g.', 'MarkerSize', 5);
plot3(PSIS_R(:,1), PSIS_R(:,2), PSIS_R(:,3), 'b.', 'MarkerSize', 5);
plot3(PSIS_L(:,1), PSIS_L(:,2), PSIS_L(:,3), 'c.', 'MarkerSize', 5);
plot3(CoM(:,1), CoM(:,2), CoM(:,3), 'k-', 'LineWidth', 2);

% Plot bounding box
vertices = [
    min_x, min_y, min_z;
    max_x, min_y, min_z;
    max_x, max_y, min_z;
    min_x, max_y, min_z;
    min_x, min_y, max_z;
    max_x, min_y, max_z;
    max_x, max_y, max_z;
    min_x, max_y, max_z;
];

faces = [
    1, 2, 3, 4;
    5, 6, 7, 8;
    1, 2, 6, 5;
    2, 3, 7, 6;
    3, 4, 8, 7;
    4, 1, 5, 8;
];

patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 1);

% Add labels and legend
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Center of Mass Trajectory');
legend('Right ASIS', 'Left ASIS', 'Right PSIS', 'Left PSIS', 'CoM', 'Bounding Box');
grid on;
axis equal;

% Plot trajectories in separate dimensions
subplot(3,2,2);
plot(1:num_frames, CoM(:,1), 'r-', 'LineWidth', 1.5);
title('CoM X Position vs. Frame');
xlabel('Frame');
ylabel('X Position (m)');
grid on;

subplot(3,2,4);
plot(1:num_frames, CoM(:,2), 'g-', 'LineWidth', 1.5);
title('CoM Y Position vs. Frame');
xlabel('Frame');
ylabel('Position (m)');
legend('Y Position');
grid on;

subplot(3,2,6)
plot(1:num_frames, CoM(:,3), 'b-', 'LineWidth', 1.5);
title('CoM Z Position vs. Frame');
xlabel('Frame');
ylabel('Position (m)');
legend('Z Position');
grid on;

%% Calculate energy cost divided in the 3 directions (X, Y, Z)

% Calculate velocities in each direction
frame_rate = 1/T.Time(2);
velocity = zeros(num_frames-1, 3);
for i = 2:num_frames
    velocity(i-1,:) = (CoM(i,:) - CoM(i-1,:)) * frame_rate;
end

% Calculate kinetic energy in each direction (1/2 * m * v^2) for each frame
kinetic_energy_x = 0.5 * mass * velocity(:,1).^2;
kinetic_energy_y = 0.5 * mass * velocity(:,2).^2;
kinetic_energy_z = 0.5 * mass * velocity(:,3).^2;
total_kinetic_energy = 0.5 * mass * sum(velocity.^2, 2);

% Calculate potential energy changes (m*g*h) for Y direction (in our case
% is the vertical!!)
gravity = 9.81; % m/s^2
potential_energy_changes_y = zeros(num_frames-1, 1);

for i = 2:num_frames
    % Y direction
    height_change_y = CoM(i,2) - CoM(i-1,2);
    potential_energy_changes_y(i-1) = mass * gravity * height_change_y;
end

% Only include positive changes in potential energy (work against gravity)
positive_potential_energy_y = potential_energy_changes_y;
positive_potential_energy_y(positive_potential_energy_y < 0) = 0;

% Calculate total mechanical work in each direction
% For X: Work = sum of absolute changes in kinetic energy
% For Y: Work = sum of absolute changes in kinetic energy + positive potential energy changes
work_x = sum(abs(diff(kinetic_energy_x)));
work_y_kinetic = sum(abs(diff(kinetic_energy_y)));
work_y_potential = sum(positive_potential_energy_y);
work_y = work_y_kinetic + work_y_potential;
work_z = sum(abs(diff(kinetic_energy_z)));

% Total work
total_work = work_x + work_y + work_z;

% Calculate percentage of energy spent in each direction
percent_work_x = (work_x / total_work) * 100;
percent_work_y = (work_y / total_work) * 100;
percent_work_z = (work_z / total_work) * 100;


%% Create figure for energy
figure('Position', [100, 100, 1000, 800]);

% chart showing work in each direction
subplot(1,2,1);
bar([work_x, work_y_kinetic, work_y_potential, work_z]);
title('Mechanical Work Distribution');
ylabel('Work (J)');
xticklabels({'X', 'Y kinetic', 'Y potential', 'Z kinetic'});
grid on;

% kinetic energy over time for each direction
subplot(1,2,2);
hold on;
plot(1:length(kinetic_energy_x), kinetic_energy_x, 'r-');
plot(1:length(kinetic_energy_y), kinetic_energy_y, 'g-');
plot(1:length(kinetic_energy_z), kinetic_energy_z, 'b-');
title('Kinetic Energy vs. Frame');
xlabel('Frame');
ylabel('Kinetic Energy (J)');
legend({'X direction', 'Y direction', 'Z direction'}, 'Location', 'best');
grid on;

%%
%% === Riepilogo Energetico ===
%% === Riepilogo Energetico ===

% Energia cinetica massima lungo Z
KE_z_max = max(kinetic_energy_z);

% Velocità massima lungo Z
v_z_max = max(abs(velocity(:,3)));

% Stampa riepilogo
fprintf('\n=== RISULTATO ENERGETICO COMPLESSIVO ===\n');
fprintf('Lavoro meccanico totale: %.2f J\n', total_work);
fprintf('  - Lavoro X (lateral):           %.2f J\n', work_x);
fprintf('  - Lavoro Y (verticale):         %.2f J (%.2f cin + %.2f pot)\n', ...
        work_y, work_y_kinetic, work_y_potential);
fprintf('  - Lavoro Z (affondo):           %.2f J\n', work_z);
fprintf('Distribuzione:\n');
fprintf('  - X: %.1f%% | Y: %.1f%% | Z: %.1f%%\n', ...
        percent_work_x, percent_work_y, percent_work_z);
fprintf('Energia cinetica massima lungo Z (affondo): %.2f J\n', KE_z_max);
fprintf('Velocità massima del CoM lungo Z: %.2f m/s\n', v_z_max);

%%
%% === Animazione nel piano ZY con picco energia Z ===
figure;
title('Traiettoria CoM nel piano ZY (con picco energia)');
xlabel('Z [m]');
ylabel('Y [m]');
axis equal;
grid on;
hold on;

% Plot statico della traiettoria
plot(CoM(:,3), CoM(:,2), 'k-', 'LineWidth', 1);

% Trova il frame del picco di energia Z
[~, idx_max_energy_Z] = max(kinetic_energy_z);

% Animazione
h_dot = plot(CoM(1,3), CoM(1,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
h_peak = plot(CoM(idx_max_energy_Z,3), CoM(idx_max_energy_Z,2), 'ro', ...
              'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Picco energia Z');

legend('Traiettoria CoM', 'CoM animato', 'Picco energia Z');

for i = 1:size(CoM,1)
    set(h_dot, 'XData', CoM(i,3), 'YData', CoM(i,2));
    pause(0.01);  % pausa tra frame
end
%%
% Trova il frame del picco di energia cinetica lungo Z
[KE_z_max, idx_max_energy_Z] = max(kinetic_energy_z);

% Tempo corrispondente (secondi)
t_max_energy_Z = T.Time(idx_max_energy_Z);

% Centro di massa in quel frame
CoM_peak = CoM(idx_max_energy_Z, :);

% Stampa risultati
fprintf('\n=== PICCO ENERGIA CINETICA LUNGO Z ===\n');
fprintf('Frame: %d\n', idx_max_energy_Z);
fprintf('Tempo: %.3f s\n', t_max_energy_Z);
fprintf('Energia cinetica massima Z: %.2f J\n', KE_z_max);
fprintf('Velocità Z nel picco: %.2f m/s\n', velocity(idx_max_energy_Z,3));
fprintf('Posizione CoM (X,Y,Z): [%.3f, %.3f, %.3f] m\n', CoM_peak(1), CoM_peak(2), CoM_peak(3));
%% confronto con le altre acquisizioni
clear; clc; close all;

% === Parametri soggetto ===
height = 1.74;     % m
mass = 73;         % kg
g = 9.81;          % m/s²
cartella = 'OneDrive_1_05-05-2025';
fileNames = {
    'Arthur trial 1.xlsx', ...
    'Arthur trial 2.xlsx', ...
    'Arthur trial 3.xlsx', ...
    'Arthur trial 4.xlsx', ...
    'Arthur trial 5.xlsx'
};

nTrials = numel(fileNames);
T_results = table('Size', [nTrials, 13], ...
    'VariableTypes', repmat({'double'}, 1, 13), ...
    'VariableNames', {'Volume_m3', 'E_cinetica_J', 'E_potenziale_J', 'E_totale_J', ...
                      'Work_X_J', 'Work_Y_J', 'Work_Z_J', ...
                      'Perc_Work_X', 'Perc_Work_Y', 'Perc_Work_Z', ...
                      'Frame_picco_Z', 'Tempo_picco_Z_s', 'Vmax_Z_m_s'});
T_results.Trial = fileNames';

for i = 1:nTrials
    fname = fileNames{i};
    fprintf('\n========= %s =========\n', fname);

    T = readtable(fullfile(cartella, fname));
    ASIS_R = [T.RightASIS_X, T.RightASIS_Y, T.RightASIS_Z];
    ASIS_L = [T.LeftASIS_X, T.LeftASIS_Y, T.LeftASIS_Z];
    PSIS_R = [T.RightPESIS_X, T.RightPESIS_Y, T.RightPESIS_Z];
    PSIS_L = [T.LeftPESIS_X, T.LeftPESIS_Y, T.LeftPESIS_Z];

    num_frames = size(ASIS_R, 1);
    CoM = zeros(num_frames, 3);
    for f = 1:num_frames
        CoM(f,:) = mean([ASIS_R(f,:); ASIS_L(f,:); PSIS_R(f,:); PSIS_L(f,:)], 1);
    end

    % === Volume occupato dal CoM ===
    range = max(CoM) - min(CoM);
    volume = prod(range);  % bounding box volume

    % === Velocità e energia ===
    frame_rate = 1 / T.Time(2);
    velocity = diff(CoM) * frame_rate;

    kinetic_energy_x = 0.5 * mass * velocity(:,1).^2;
    kinetic_energy_y = 0.5 * mass * velocity(:,2).^2;
    kinetic_energy_z = 0.5 * mass * velocity(:,3).^2;
    total_kinetic_energy = kinetic_energy_x + kinetic_energy_y + kinetic_energy_z;

    % === Energia potenziale (solo nei tratti in salita) ===
    delta_h = diff(CoM(:,2));
    potential_energy_changes = mass * g * delta_h;
    potential_energy_changes(potential_energy_changes < 0) = 0;
    potential_energy = sum(potential_energy_changes);

    % === Lavoro meccanico ===
    work_x = sum(abs(diff(kinetic_energy_x)));
    work_y_kinetic = sum(abs(diff(kinetic_energy_y)));
    work_y = work_y_kinetic + potential_energy;
    work_z = sum(abs(diff(kinetic_energy_z)));
    total_work = work_x + work_y + work_z;

    % === Percentuali ===
    perc_x = (work_x / total_work) * 100;
    perc_y = (work_y / total_work) * 100;
    perc_z = (work_z / total_work) * 100;

    % === Picco energia Z ===
    [E_k_max_Z, idx_max_Z] = max(kinetic_energy_z);
    t_max = T.Time(idx_max_Z);
    v_z_max = velocity(idx_max_Z,3);

    % === Salvataggio risultati ===
    T_results.Volume_m3(i)       = volume;
    T_results.E_cinetica_J(i)    = sum(total_kinetic_energy);
    T_results.E_potenziale_J(i)  = potential_energy;
    T_results.E_totale_J(i)      = sum(total_kinetic_energy) + potential_energy;
    T_results.Work_X_J(i)        = work_x;
    T_results.Work_Y_J(i)        = work_y;
    T_results.Work_Z_J(i)        = work_z;
    T_results.Perc_Work_X(i)     = perc_x;
    T_results.Perc_Work_Y(i)     = perc_y;
    T_results.Perc_Work_Z(i)     = perc_z;
    T_results.Frame_picco_Z(i)   = idx_max_Z;
    T_results.Tempo_picco_Z_s(i) = t_max;
    T_results.Vmax_Z_m_s(i)      = v_z_max;

    % === Output console ===
    fprintf('Volume: %.6f m³ | Energia meccanica: %.2f J\n', volume, T_results.E_totale_J(i));
    fprintf('Work → X: %.2f J, Y: %.2f J, Z: %.2f J\n', work_x, work_y, work_z);
    fprintf('Distribuzione: X %.1f%% | Y %.1f%% | Z %.1f%%\n', perc_x, perc_y, perc_z);
    fprintf('Picco Ek Z @ frame %d (t = %.3f s) | Vz max = %.2f m/s\n', idx_max_Z, t_max, v_z_max);

    % === Animazione ZY ===
    figure('Name', fname);
    hold on; grid on; axis equal;
    title(['CoM ZY - ', fname]);
    xlabel('Z [m]'); ylabel('Y [m]');
    plot(CoM(:,3), CoM(:,2), 'k-', 'LineWidth', 1);
    plot(CoM(idx_max_Z,3), CoM(idx_max_Z,2), 'ro', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r');

    h_dot = plot(CoM(1,3), CoM(1,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    for f = 1:num_frames
        set(h_dot, 'XData', CoM(f,3), 'YData', CoM(f,2));
        pause(0.01);
    end
end

% === Tabella finale ===
disp('✅ RIEPILOGO COMPLETO:');
disp(T_results);
