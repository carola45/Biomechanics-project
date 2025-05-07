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

    %% Volume e traiettoria del CoM
    % Calcolare il volume occupato dal CoM (ginocchio)
    range = max(CoM_knee) - min(CoM_knee);
    volume_knee = prod(range);  % bounding box volume

    % === Visualizzazione 3D del volume del CoM del ginocchio ===
    figure('Name', ['3D Volume CoM Ginocchio - ', fname]);
    hold on; grid on; axis equal;
    title(['Volume occupato dal CoM del Ginocchio - ', fname]);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');

    % Calcola i vertici del bounding box
    min_coords_knee = min(CoM_knee);
    max_coords_knee = max(CoM_knee);
    
    % Vertici del cubo che rappresenta il volume
    vertices_knee = [
        min_coords_knee(1), min_coords_knee(2), min_coords_knee(3);
        max_coords_knee(1), min_coords_knee(2), min_coords_knee(3);
        max_coords_knee(1), max_coords_knee(2), min_coords_knee(3);
        min_coords_knee(1), max_coords_knee(2), min_coords_knee(3);
        min_coords_knee(1), min_coords_knee(2), max_coords_knee(3);
        max_coords_knee(1), min_coords_knee(2), max_coords_knee(3);
        max_coords_knee(1), max_coords_knee(2), max_coords_knee(3);
        min_coords_knee(1), max_coords_knee(2), max_coords_knee(3)
    ];

    % Facce del cubo (triangoli per il plot)
    faces_knee = [
        1 2 3 4;
        5 6 7 8;
        1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        4 1 5 8
    ];

    % Disegna il volume del ginocchio
    patch('Vertices', vertices_knee, 'Faces', faces_knee, 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'k');

    % Plotta la traiettoria del ginocchio (CoM)
    plot3(CoM_knee(:,1), CoM_knee(:,2), CoM_knee(:,3), 'k-', 'LineWidth', 1.5);
    scatter3(CoM_knee(:,1), CoM_knee(:,2), CoM_knee(:,3), 30, 'r', 'filled');

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
end
