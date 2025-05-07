clear all
clc

% Parametri soggetto
height = 1.74;       % Altezza (m)
mass = 73;           % Massa (kg)

filename = 'Arthur trial 2.xlsx';
T = readtable(filename);

% Estrazione delle posizioni dei marker del ginocchio
lateral_epicondyle = [T{:,36}, T{:,37}, T{:,38}];  % Right lateral epicondyle of femur (X, Y, Z)
fibula_head = [T{:,39}, T{:,40}, T{:,41}];          % Right head of the fibula (X, Y, Z)
medial_epicondyle = [T{:,41}, T{:,42}, T{:,43}];    % Right medial epicondyle of femur (X, Y, Z)

% Numero di frame
num_frames = size(lateral_epicondyle, 1);

% Calcolare il centroide (CoM) del triangolo formato dai marker per ogni frame
CoM_knee = zeros(num_frames, 3);
for f = 1:num_frames
    CoM_knee(f,:) = (lateral_epicondyle(f,:) + fibula_head(f,:) + medial_epicondyle(f,:)) / 3;
end

% === Volume occupato dal CoM ===
range = max(CoM_knee) - min(CoM_knee);
volume = prod(range);  % bounding box volume

% === Visualizzazione 3D con X rappresentato lungo Z e Z lungo X ===
figure('Position', [100, 100, 1200, 800]);

% Plot dei marker e della traiettoria del CoM del ginocchio
subplot(3,2,[1,3,5]);
hold on;
grid on;
axis equal;

% Cambia le coordinate per gli assi come richiesto
% Assegniamo le coordinate in modo che l'asse X sia lungo Z e l'asse Z lungo X
plot3(CoM_knee(:,3), CoM_knee(:,2), CoM_knee(:,1), 'k-', 'LineWidth', 2);
scatter3(CoM_knee(:,3), CoM_knee(:,2), CoM_knee(:,1), 30, 'r', 'filled');

% Legenda e etichette
xlabel('Z (m)'); ylabel('Y (m)'); zlabel('X (m)');
title('Traiettoria CoM del Ginocchio (assi modificati)');
legend('Traiettoria CoM', 'Punti del CoM');
grid on;

% === Animazione della traiettoria in 3D ===
h_dot = plot3(CoM_knee(1,3), CoM_knee(1,2), CoM_knee(1,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
h_path = plot3([], [], [], 'r-', 'LineWidth', 1.5);

for f = 1:num_frames
    % Aggiorna la posizione del punto del CoM
    set(h_dot, 'XData', CoM_knee(f,3), 'YData', CoM_knee(f,2), 'ZData', CoM_knee(f,1));

    % Aggiorna la traiettoria
    x_data = CoM_knee(1:f, 3);
    y_data = CoM_knee(1:f, 2);
    z_data = CoM_knee(1:f, 1);
    set(h_path, 'XData', x_data, 'YData', y_data, 'ZData', z_data);

    % Pausa per l'animazione
    pause(0.01);
end

% === Output del volume del CoM ===
disp('Volume occupato dal CoM del Ginocchio:');
disp(volume);
%%
clear all
clc

% Parametri soggetto
height = 1.74;       % Altezza (m)
mass = 73;           % Massa (kg)

filename = 'Arthur trial 2.xlsx';
T = readtable(filename);

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

% Calcolare la distanza percorsa solo fino al massimo di flessione (senza ritorno)
[~, idx_max_Z] = min(CoM_knee(:,3));  % Trova il punto di massima flessione (minimo lungo Z)

% Calcolare la distanza percorsa fino al punto di massima flessione
total_distance = 0;
for f = 2:idx_max_Z
    % Calcolare la distanza tra il frame attuale e il precedente
    segment_distance = norm(CoM_knee(f,:) - CoM_knee(f-1,:));
    total_distance = total_distance + segment_distance;
end

% Output della distanza percorsa
fprintf('La distanza percorsa dal ginocchio fino alla massima flessione (escludendo il ritorno) è: %.3f m\n', total_distance);

%% Animazione 3D con visualizzazione della traiettoria e dei punti di interesse

figure;
hold on;
grid on;
axis equal;
title('Traiettoria CoM del ginocchio fino alla massima flessione');
xlabel('Z [m]');
ylabel('Y [m]');
zlabel('X [m]');

% Plot della traiettoria del ginocchio (CoM)
plot3(CoM_knee(:,1), CoM_knee(:,2), CoM_knee(:,3), 'k-', 'LineWidth', 1.5);

% Aggiungi il punto iniziale
plot3(CoM_knee(1,1), CoM_knee(1,2), CoM_knee(1,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Posizione Iniziale');

% Aggiungi il punto di massimo Z (massima flessione)
plot3(CoM_knee(idx_max_Z, 1), CoM_knee(idx_max_Z, 2), CoM_knee(idx_max_Z, 3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Massima Flessione');

% Animazione fino al punto di massima flessione
h_dot = plot3(CoM_knee(1,1), CoM_knee(1,2), CoM_knee(1,3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
for f = 1:idx_max_Z
    set(h_dot, 'XData', CoM_knee(f,1), 'YData', CoM_knee(f,2), 'ZData', CoM_knee(f,3));
    pause(0.01);
end

% Plot della distanza percorsa (linea tra il punto iniziale e il massimo Z)
plot3([CoM_knee(1,1), CoM_knee(idx_max_Z,1)], [CoM_knee(1,2), CoM_knee(idx_max_Z,2)], [CoM_knee(1,3), CoM_knee(idx_max_Z,3)], 'm-', 'LineWidth', 2, 'DisplayName', 'Distanza Percorsa');

legend;


