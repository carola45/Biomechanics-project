% S1_recontruction
% 1-2-3)  missed markers (half of data): 
%	- Left medial epicondyle of femur (X,Y,Z) 
%
% 5)	missed markers (half of data):
% 	- Right Thigh


%% trial 1

filename = 'SubjectBis 1 trial 1.xlsx';
T = readtable(filename);

% Extract marker positions
P1 = [T.LeftTrochanter_X, T.LeftTrochanter_Y, T.LeftTrochanter_Z];
P2 = [T.LeftLateralEpicondyleOfFemur_X, T.LeftLateralEpicondyleOfFemur_Y, T.LeftLateralEpicondyleOfFemur_Z];
P3 = [T.LeftThigh_X, T.LeftThigh_Y, T.LeftThigh_Z];
P4 = [T.LeftMedialEpicondyleOfFemur_X, T.LeftMedialEpicondyleOfFemur_Y, T.LeftMedialEpicondyleOfFemur_Z];

% Find frames where P4 is valid and where it's missing
valid_idx = ~any(isnan(P4), 2);
missing_idx = find(~valid_idx);

% Create copy for reconstruction
P4_recon = P4;

% use local coordinate system that evolves with movement
% For each missing frame, find nearest valid frames (before and after)
for i = missing_idx'
    % Find nearest valid frames before and after current frame
    prev_valid = find(valid_idx & (1:length(valid_idx))' < i, 1, 'last');
    next_valid = find(valid_idx & (1:length(valid_idx))' > i, 1, 'first');
    
    % If we have both bounds, interpolate between them
    if ~isempty(prev_valid) && ~isempty(next_valid)
        % Calculate weights based on temporal distance
        total_dist = next_valid - prev_valid;
        weight_next = (i - prev_valid) / total_dist;
        weight_prev = 1 - weight_next;
        
        % Get local coordinates in both reference frames
        local_prev = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        local_next = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Interpolate local coordinates
        local_interp = weight_prev * local_prev + weight_next * local_next;
        
        % Convert back to global coordinates using current frame's reference markers
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_interp);
    
    % If we only have previous valid frames
    elseif ~isempty(prev_valid)
        % Get local coordinates in previous valid frame
        local_coords = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If we only have future valid frames
    elseif ~isempty(next_valid)
        % Get local coordinates in next valid frame
        local_coords = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If no valid frames available (unlikely but possible)
    else
        % Fall back to the sphere intersection method
        % Calculate mean distances from valid frames
        d1 = mean(sqrt(sum((P4(valid_idx,:) - P1(valid_idx,:)).^2, 2)));
        d2 = mean(sqrt(sum((P4(valid_idx,:) - P2(valid_idx,:)).^2, 2)));
        d3 = mean(sqrt(sum((P4(valid_idx,:) - P3(valid_idx,:)).^2, 2)));
        
        % Known markers positions at frame i
        p1 = P1(i,:);
        p2 = P2(i,:);
        p3 = P3(i,:);
        
        % Create local coordinate system
        ex = (p2 - p1) / norm(p2 - p1);
        i_val = dot(ex, p3 - p1);
        ey = (p3 - p1 - i_val*ex);
        ey = ey / norm(ey);
        ez = cross(ex, ey);
        d = norm(p2 - p1);
        j_val = dot(ey, p3 - p1);
        
        % Calculate coordinates in local system
        x = (d1^2 - d2^2 + d^2) / (2*d);
        y = (d1^2 - d3^2 + i_val^2 + j_val^2 - 2*i_val*x) / (2*j_val);
        z_sq = d1^2 - x^2 - y^2;
        
        % Ensure we don't try to take square root of negative number
        if z_sq >= 0
            z = sqrt(z_sq);
            P4_recon(i,:) = p1 + x*ex + y*ey + z*ez;
        else
            % Handle case where there's numerical error
            z = 0;
            P4_recon(i,:) = p1 + x*ex + y*ey;
        end
    end
end

% Apply a small amount of smoothing to ensure continuity at transitions
window_size = 5; % Adjust based on sampling rate and movement dynamics
P4_recon = smoothTrajectory(P4_recon, valid_idx, window_size);

% Plot results
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');


% Plot known trajectories
plot3(P1(:,1), P1(:,2), P1(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'Left Trochanter');
plot3(P2(:,1), P2(:,2), P2(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Left Lateral Epicondyle');
plot3(P3(:,1), P3(:,2), P3(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Left Thigh');

% Plot original P4 on valid frames
plot3(P4(valid_idx,1), P4(valid_idx,2), P4(valid_idx,3), 'ko', 'MarkerSize', 4, 'DisplayName', 'Original P4 (Valid)');

% Plot reconstructed P4 on missing frames
plot3(P4_recon(missing_idx,1), P4_recon(missing_idx,2), P4_recon(missing_idx,3), 'mx', 'MarkerSize', 6, 'DisplayName', 'Reconstructed P4');

% Plot the full reconstructed trajectory
plot3(P4_recon(:,1), P4_recon(:,2), P4_recon(:,3), 'm--', 'LineWidth', 1, 'DisplayName', 'Full Reconstructed Path');

legend('Location', 'best');

% Create a new Excel file with gaps filled
% Make a copy of the original table
T_filled = T;

% Replace the missing values with reconstructed ones
T_filled.LeftMedialEpicondyleOfFemur_X(missing_idx) = P4_recon(missing_idx, 1);
T_filled.LeftMedialEpicondyleOfFemur_Y(missing_idx) = P4_recon(missing_idx, 2);
T_filled.LeftMedialEpicondyleOfFemur_Z(missing_idx) = P4_recon(missing_idx, 3);

% Create output filename
[filepath, name, ext] = fileparts(filename);
output_filename = fullfile(filepath, [name '_filled' ext]);

% Write the filled data to a new Excel file
writetable(T_filled, output_filename);
fprintf('Filled data saved to: %s\n', output_filename);


% Helper function to get local coordinates
function local_coords = getLocalCoordinates(p1, p2, p3, p4)
    % Create local coordinate system
    ex = (p2 - p1) / norm(p2 - p1);
    temp = p3 - p1;
    i_val = dot(ex, temp);
    ey = temp - i_val*ex;
    ey = ey / norm(ey);
    ez = cross(ex, ey);
    
    % Transform p4 into local coordinates
    vec = p4 - p1;
    local_coords = [dot(vec, ex), dot(vec, ey), dot(vec, ez)];
end

% Helper function to transform back to global coordinates
function global_coords = transformToGlobal(p1, p2, p3, local_coords)
    % Create local coordinate system
    ex = (p2 - p1) / norm(p2 - p1);
    temp = p3 - p1;
    i_val = dot(ex, temp);
    ey = temp - i_val*ex;
    ey = ey / norm(ey);
    ez = cross(ex, ey);
    
    % Transform from local to global
    global_coords = p1 + local_coords(1)*ex + local_coords(2)*ey + local_coords(3)*ez;
end

% Helper function to smooth trajectory across transitions
function smoothed = smoothTrajectory(trajectory, valid_idx, window_size)
    smoothed = trajectory;
    
    % Only smooth around transitions between original and reconstructed data
    transitions = find(diff([0; valid_idx; 0]) ~= 0);
    
    for i = 1:length(transitions)
        center = transitions(i);
        % Define window boundaries, respecting array limits
        start_idx = max(1, center - floor(window_size/2));
        end_idx = min(size(trajectory, 1), center + floor(window_size/2));
        
        % Apply smoothing only if we have enough points
        if end_idx - start_idx >= 3
            for dim = 1:3
                % Use moving average smoothing
                smoothed(start_idx:end_idx, dim) = movmean(trajectory(start_idx:end_idx, dim), 3, 'omitnan');
            end
        end
    end
end

%% trial 2

filename = 'SubjectBis 1 trial 2.xlsx';
T = readtable(filename);

% Extract marker positions
P1 = [T.LeftTrochanter_X, T.LeftTrochanter_Y, T.LeftTrochanter_Z];
P2 = [T.LeftLateralEpicondyleOfFemur_X, T.LeftLateralEpicondyleOfFemur_Y, T.LeftLateralEpicondyleOfFemur_Z];
P3 = [T.LeftThigh_X, T.LeftThigh_Y, T.LeftThigh_Z];
P4 = [T.LeftMedialEpicondyleOfFemur_X, T.LeftMedialEpicondyleOfFemur_Y, T.LeftMedialEpicondyleOfFemur_Z];

% Find frames where P4 is valid and where it's missing
valid_idx = ~any(isnan(P4), 2);
missing_idx = find(~valid_idx);

% Create copy for reconstruction
P4_recon = P4;

% use local coordinate system that evolves with movement
% For each missing frame, find nearest valid frames (before and after)
for i = missing_idx'
    % Find nearest valid frames before and after current frame
    prev_valid = find(valid_idx & (1:length(valid_idx))' < i, 1, 'last');
    next_valid = find(valid_idx & (1:length(valid_idx))' > i, 1, 'first');
    
    % If we have both bounds, interpolate between them
    if ~isempty(prev_valid) && ~isempty(next_valid)
        % Calculate weights based on temporal distance
        total_dist = next_valid - prev_valid;
        weight_next = (i - prev_valid) / total_dist;
        weight_prev = 1 - weight_next;
        
        % Get local coordinates in both reference frames
        local_prev = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        local_next = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Interpolate local coordinates
        local_interp = weight_prev * local_prev + weight_next * local_next;
        
        % Convert back to global coordinates using current frame's reference markers
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_interp);
    
    % If we only have previous valid frames
    elseif ~isempty(prev_valid)
        % Get local coordinates in previous valid frame
        local_coords = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If we only have future valid frames
    elseif ~isempty(next_valid)
        % Get local coordinates in next valid frame
        local_coords = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If no valid frames available (unlikely but possible)
    else
        % Fall back to the sphere intersection method
        % Calculate mean distances from valid frames
        d1 = mean(sqrt(sum((P4(valid_idx,:) - P1(valid_idx,:)).^2, 2)));
        d2 = mean(sqrt(sum((P4(valid_idx,:) - P2(valid_idx,:)).^2, 2)));
        d3 = mean(sqrt(sum((P4(valid_idx,:) - P3(valid_idx,:)).^2, 2)));
        
        % Known markers positions at frame i
        p1 = P1(i,:);
        p2 = P2(i,:);
        p3 = P3(i,:);
        
        % Create local coordinate system
        ex = (p2 - p1) / norm(p2 - p1);
        i_val = dot(ex, p3 - p1);
        ey = (p3 - p1 - i_val*ex);
        ey = ey / norm(ey);
        ez = cross(ex, ey);
        d = norm(p2 - p1);
        j_val = dot(ey, p3 - p1);
        
        % Calculate coordinates in local system
        x = (d1^2 - d2^2 + d^2) / (2*d);
        y = (d1^2 - d3^2 + i_val^2 + j_val^2 - 2*i_val*x) / (2*j_val);
        z_sq = d1^2 - x^2 - y^2;
        
        % Ensure we don't try to take square root of negative number
        if z_sq >= 0
            z = sqrt(z_sq);
            P4_recon(i,:) = p1 + x*ex + y*ey + z*ez;
        else
            % Handle case where there's numerical error
            z = 0;
            P4_recon(i,:) = p1 + x*ex + y*ey;
        end
    end
end

% Apply a small amount of smoothing to ensure continuity at transitions
window_size = 5; % Adjust based on sampling rate and movement dynamics
P4_recon = smoothTrajectory(P4_recon, valid_idx, window_size);

% Plot results
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

% Plot known trajectories
plot3(P1(:,1), P1(:,2), P1(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'Left Trochanter');
plot3(P2(:,1), P2(:,2), P2(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Left Lateral Epicondyle');
plot3(P3(:,1), P3(:,2), P3(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Left Thigh');

% Plot original P4 on valid frames
plot3(P4(valid_idx,1), P4(valid_idx,2), P4(valid_idx,3), 'ko', 'MarkerSize', 4, 'DisplayName', 'Original P4 (Valid)');

% Plot reconstructed P4 on missing frames
plot3(P4_recon(missing_idx,1), P4_recon(missing_idx,2), P4_recon(missing_idx,3), 'mx', 'MarkerSize', 6, 'DisplayName', 'Reconstructed P4');

% Plot the full reconstructed trajectory
plot3(P4_recon(:,1), P4_recon(:,2), P4_recon(:,3), 'm--', 'LineWidth', 1, 'DisplayName', 'Full Reconstructed Path');

legend('Location', 'best');

% Create a new Excel file with gaps filled
% Make a copy of the original table
T_filled = T;

% Replace the missing values with reconstructed ones
T_filled.LeftMedialEpicondyleOfFemur_X(missing_idx) = P4_recon(missing_idx, 1);
T_filled.LeftMedialEpicondyleOfFemur_Y(missing_idx) = P4_recon(missing_idx, 2);
T_filled.LeftMedialEpicondyleOfFemur_Z(missing_idx) = P4_recon(missing_idx, 3);

% Create output filename
[filepath, name, ext] = fileparts(filename);
output_filename = fullfile(filepath, [name '_filled' ext]);

% Write the filled data to a new Excel file
writetable(T_filled, output_filename);
fprintf('Filled data saved to: %s\n', output_filename);

%% trial 3

filename = 'SubjectBis 1 trial 3.xlsx';
T = readtable(filename);

% Extract marker positions
P1 = [T.LeftTrochanter_X, T.LeftTrochanter_Y, T.LeftTrochanter_Z];
P2 = [T.LeftLateralEpicondyleOfFemur_X, T.LeftLateralEpicondyleOfFemur_Y, T.LeftLateralEpicondyleOfFemur_Z];
P3 = [T.LeftThigh_X, T.LeftThigh_Y, T.LeftThigh_Z];
P4 = [T.LeftMedialEpicondyleOfFemur_X, T.LeftMedialEpicondyleOfFemur_Y, T.LeftMedialEpicondyleOfFemur_Z];

% Find frames where P4 is valid and where it's missing
valid_idx = ~any(isnan(P4), 2);
missing_idx = find(~valid_idx);

% Create copy for reconstruction
P4_recon = P4;

% use local coordinate system that evolves with movement
% For each missing frame, find nearest valid frames (before and after)
for i = missing_idx'
    % Find nearest valid frames before and after current frame
    prev_valid = find(valid_idx & (1:length(valid_idx))' < i, 1, 'last');
    next_valid = find(valid_idx & (1:length(valid_idx))' > i, 1, 'first');
    
    % If we have both bounds, interpolate between them
    if ~isempty(prev_valid) && ~isempty(next_valid)
        % Calculate weights based on temporal distance
        total_dist = next_valid - prev_valid;
        weight_next = (i - prev_valid) / total_dist;
        weight_prev = 1 - weight_next;
        
        % Get local coordinates in both reference frames
        local_prev = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        local_next = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Interpolate local coordinates
        local_interp = weight_prev * local_prev + weight_next * local_next;
        
        % Convert back to global coordinates using current frame's reference markers
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_interp);
    
    % If we only have previous valid frames
    elseif ~isempty(prev_valid)
        % Get local coordinates in previous valid frame
        local_coords = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If we only have future valid frames
    elseif ~isempty(next_valid)
        % Get local coordinates in next valid frame
        local_coords = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If no valid frames available (unlikely but possible)
    else
        % Fall back to the sphere intersection method
        % Calculate mean distances from valid frames
        d1 = mean(sqrt(sum((P4(valid_idx,:) - P1(valid_idx,:)).^2, 2)));
        d2 = mean(sqrt(sum((P4(valid_idx,:) - P2(valid_idx,:)).^2, 2)));
        d3 = mean(sqrt(sum((P4(valid_idx,:) - P3(valid_idx,:)).^2, 2)));
        
        % Known markers positions at frame i
        p1 = P1(i,:);
        p2 = P2(i,:);
        p3 = P3(i,:);
        
        % Create local coordinate system
        ex = (p2 - p1) / norm(p2 - p1);
        i_val = dot(ex, p3 - p1);
        ey = (p3 - p1 - i_val*ex);
        ey = ey / norm(ey);
        ez = cross(ex, ey);
        d = norm(p2 - p1);
        j_val = dot(ey, p3 - p1);
        
        % Calculate coordinates in local system
        x = (d1^2 - d2^2 + d^2) / (2*d);
        y = (d1^2 - d3^2 + i_val^2 + j_val^2 - 2*i_val*x) / (2*j_val);
        z_sq = d1^2 - x^2 - y^2;
        
        % Ensure we don't try to take square root of negative number
        if z_sq >= 0
            z = sqrt(z_sq);
            P4_recon(i,:) = p1 + x*ex + y*ey + z*ez;
        else
            % Handle case where there's numerical error
            z = 0;
            P4_recon(i,:) = p1 + x*ex + y*ey;
        end
    end
end

% Apply a small amount of smoothing to ensure continuity at transitions
window_size = 5; % Adjust based on sampling rate and movement dynamics
P4_recon = smoothTrajectory(P4_recon, valid_idx, window_size);

% Plot results
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

% Plot known trajectories
plot3(P1(:,1), P1(:,2), P1(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'Left Trochanter');
plot3(P2(:,1), P2(:,2), P2(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Left Lateral Epicondyle');
plot3(P3(:,1), P3(:,2), P3(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Left Thigh');

% Plot original P4 on valid frames
plot3(P4(valid_idx,1), P4(valid_idx,2), P4(valid_idx,3), 'ko', 'MarkerSize', 4, 'DisplayName', 'Original P4 (Valid)');

% Plot reconstructed P4 on missing frames
plot3(P4_recon(missing_idx,1), P4_recon(missing_idx,2), P4_recon(missing_idx,3), 'mx', 'MarkerSize', 6, 'DisplayName', 'Reconstructed P4');

% Plot the full reconstructed trajectory
plot3(P4_recon(:,1), P4_recon(:,2), P4_recon(:,3), 'm--', 'LineWidth', 1, 'DisplayName', 'Full Reconstructed Path');

legend('Location', 'best');

% Create a new Excel file with gaps filled
% Make a copy of the original table
T_filled = T;

% Replace the missing values with reconstructed ones
T_filled.LeftMedialEpicondyleOfFemur_X(missing_idx) = P4_recon(missing_idx, 1);
T_filled.LeftMedialEpicondyleOfFemur_Y(missing_idx) = P4_recon(missing_idx, 2);
T_filled.LeftMedialEpicondyleOfFemur_Z(missing_idx) = P4_recon(missing_idx, 3);

% Create output filename
[filepath, name, ext] = fileparts(filename);
output_filename = fullfile(filepath, [name '_filled' ext]);

% Write the filled data to a new Excel file
writetable(T_filled, output_filename);
fprintf('Filled data saved to: %s\n', output_filename);



%% trial 5 - same approach

filename = 'SubjectBis 1 trial 5.xlsx';
T = readtable(filename);
P1 = [T.RightTrochanter_X, T.RightTrochanter_Y, T.RightTrochanter_Z];
P2 = [T.RightLateralEpicondyleOfFemur_X, T.RightLateralEpicondyleOfFemur_Y, T.RightLateralEpicondyleOfFemur_Z];
P3 = [T.RightMedialEpicondyleOfFemur_X, T.RightMedialEpicondyleOfFemur_Y, T.RightMedialEpicondyleOfFemur_Z];

P4 = [T.RightThigh_X, T.RightThigh_Y, T.RightThigh_Z];

% Find frames where P4 is valid and where it's missing
valid_idx = ~any(isnan(P4), 2);
missing_idx = find(~valid_idx);

% Create copy for reconstruction
P4_recon = P4;

% use local coordinate system that evolves with movement
% For each missing frame, find nearest valid frames (before and after)
for i = missing_idx'
    % Find nearest valid frames before and after current frame
    prev_valid = find(valid_idx & (1:length(valid_idx))' < i, 1, 'last');
    next_valid = find(valid_idx & (1:length(valid_idx))' > i, 1, 'first');
    
    % If we have both bounds, interpolate between them
    if ~isempty(prev_valid) && ~isempty(next_valid)
        % Calculate weights based on temporal distance
        total_dist = next_valid - prev_valid;
        weight_next = (i - prev_valid) / total_dist;
        weight_prev = 1 - weight_next;
        
        % Get local coordinates in both reference frames
        local_prev = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        local_next = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Interpolate local coordinates
        local_interp = weight_prev * local_prev + weight_next * local_next;
        
        % Convert back to global coordinates using current frame's reference markers
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_interp);
    
    % If we only have previous valid frames
    elseif ~isempty(prev_valid)
        % Get local coordinates in previous valid frame
        local_coords = getLocalCoordinates(P1(prev_valid,:), P2(prev_valid,:), P3(prev_valid,:), P4(prev_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If we only have future valid frames
    elseif ~isempty(next_valid)
        % Get local coordinates in next valid frame
        local_coords = getLocalCoordinates(P1(next_valid,:), P2(next_valid,:), P3(next_valid,:), P4(next_valid,:));
        
        % Apply to current frame
        P4_recon(i,:) = transformToGlobal(P1(i,:), P2(i,:), P3(i,:), local_coords);
    
    % If no valid frames available (unlikely but possible)
    else
        % Fall back to the sphere intersection method
        % Calculate mean distances from valid frames
        d1 = mean(sqrt(sum((P4(valid_idx,:) - P1(valid_idx,:)).^2, 2)));
        d2 = mean(sqrt(sum((P4(valid_idx,:) - P2(valid_idx,:)).^2, 2)));
        d3 = mean(sqrt(sum((P4(valid_idx,:) - P3(valid_idx,:)).^2, 2)));
        
        % Known markers positions at frame i
        p1 = P1(i,:);
        p2 = P2(i,:);
        p3 = P3(i,:);
        
        % Create local coordinate system
        ex = (p2 - p1) / norm(p2 - p1);
        i_val = dot(ex, p3 - p1);
        ey = (p3 - p1 - i_val*ex);
        ey = ey / norm(ey);
        ez = cross(ex, ey);
        d = norm(p2 - p1);
        j_val = dot(ey, p3 - p1);
        
        % Calculate coordinates in local system
        x = (d1^2 - d2^2 + d^2) / (2*d);
        y = (d1^2 - d3^2 + i_val^2 + j_val^2 - 2*i_val*x) / (2*j_val);
        z_sq = d1^2 - x^2 - y^2;
        
        % Ensure we don't try to take square root of negative number
        if z_sq >= 0
            z = sqrt(z_sq);
            P4_recon(i,:) = p1 + x*ex + y*ey + z*ez;
        else
            % Handle case where there's numerical error
            z = 0;
            P4_recon(i,:) = p1 + x*ex + y*ey;
        end
    end
end

% Apply a small amount of smoothing to ensure continuity at transitions
window_size = 5; % Adjust based on sampling rate and movement dynamics
P4_recon = smoothTrajectory(P4_recon, valid_idx, window_size);

% Plot results
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

% Plot known trajectories
plot3(P1(:,1), P1(:,2), P1(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'Left Trochanter');
plot3(P2(:,1), P2(:,2), P2(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Left Lateral Epicondyle');
plot3(P3(:,1), P3(:,2), P3(:,3), 'b-', 'LineWidth', 1, 'DisplayName', 'Left Thigh');

% Plot original P4 on valid frames
plot3(P4(valid_idx,1), P4(valid_idx,2), P4(valid_idx,3), 'ko', 'MarkerSize', 4, 'DisplayName', 'Original P4 (Valid)');

% Plot reconstructed P4 on missing frames
plot3(P4_recon(missing_idx,1), P4_recon(missing_idx,2), P4_recon(missing_idx,3), 'mx', 'MarkerSize', 6, 'DisplayName', 'Reconstructed P4');

% Plot the full reconstructed trajectory
plot3(P4_recon(:,1), P4_recon(:,2), P4_recon(:,3), 'm--', 'LineWidth', 1, 'DisplayName', 'Full Reconstructed Path');

legend('Location', 'best');

% Create a new Excel file with gaps filled
% Make a copy of the original table
T_filled = T;

% Replace the missing values with reconstructed ones
T_filled.RightThigh_X(missing_idx) = P4_recon(missing_idx, 1);
T_filled.RightThigh_Y(missing_idx) = P4_recon(missing_idx, 2);
T_filled.RightThigh_Z(missing_idx) = P4_recon(missing_idx, 3);

% Create output filename
[filepath, name, ext] = fileparts(filename);
output_filename = fullfile(filepath, [name '_filled' ext]);

% Write the filled data to a new Excel file
writetable(T_filled, output_filename);
fprintf('Filled data saved to: %s\n', output_filename);

