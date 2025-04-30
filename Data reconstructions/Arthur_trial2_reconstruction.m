% recostruction of data
% arthur prova 2 : LEFT THIGH marker is missing

% I tried to recontruct the missing positions by imposing the distance of the
% missing marker from the other markers on the same limb
% the imposed distance is the mean of the distance of the markers in the
% frames were the data of the missing marker is avaiable. Feel free to
% correct

filename = 'Arthur trial 2.xlsx';
T = readtable(filename);
P1 = [T.LeftTrochanter_X, T.LeftTrochanter_Y, T.LeftTrochanter_Z];
P2 = [T.LeftLateralEpicondyleOfFemur_X, T.LeftLateralEpicondyleOfFemur_Y, T.LeftLateralEpicondyleOfFemur_Z];
P3 = [T.LeftMedialEpicondyleOfFemur_X, T.LeftMedialEpicondyleOfFemur_Y, T.LeftMedialEpicondyleOfFemur_Z];

P4 = [T.LeftThigh_X, T.LeftThigh_Y, T.LeftThigh_Z];

% find frame were P4 is valid
valid_idx = ~any(isnan(P4), 2);
missing_idx = ~valid_idx;

% mean distances from P4 to each marker over valid frames
d1 = mean(sqrt(sum((P4(valid_idx,:) - P1(valid_idx,:)).^2, 2)));
d2 = mean(sqrt(sum((P4(valid_idx,:) - P2(valid_idx,:)).^2, 2)));
d3 = mean(sqrt(sum((P4(valid_idx,:) - P3(valid_idx,:)).^2, 2)));

% initialize reconstructed P4
P4_recon = P4;

% reconstruction for missing frames
for i = find(missing_idx)'
    % known markers positions at frame i
    p1 = P1(i,:);
    p2 = P2(i,:);
    p3 = P3(i,:);

    % intersection of three sferes centered in the 3 known markers, with
    % radius equal to the mean distance
    % 
    % definisci il sistema locale
    ex = (p2 - p1) / norm(p2 - p1);
    i_val = dot(ex, p3 - p1);
    ey = (p3 - p1 - i_val*ex);
    ey = ey / norm(ey);
    ez = cross(ex, ey);
    
    d = norm(p2 - p1);
    j_val = dot(ey, p3 - p1);
    
    % coordinate nel sistema locale
    x = (d1^2 - d2^2 + d^2) / (2*d);
    y = (d1^2 - d3^2 + i_val^2 + j_val^2 - 2*i_val*x) / (2*j_val);
    
    z_sq = d1^2 - x^2 - y^2;
    
    z = -sqrt(z_sq);  % una delle due soluzioni
    
    % costruisci punto nel sistema globale
    P = p1 + x*ex + y*ey + z*ez;

    P4_recon(i,:) = P;
end

%% plot
figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Reconstruction of P4');

% Plot known trajectories
plot3(P1(:,1),P1(:,2),P1(:,3), 'r-');
plot3(P2(:,1),P2(:,2),P2(:,3), 'g-');
plot3(P3(:,1),P3(:,2),P3(:,3), 'b-');

% Plot original P4 on valid frames
plot3(P4(valid_idx,1), P4(valid_idx,2), P4(valid_idx,3), 'ko');
% Plot reconstructed P4 on missing frames
missing_idx = ~valid_idx;
plot3(P4_recon(missing_idx,1), P4_recon(missing_idx,2), P4_recon(missing_idx,3), 'mx');





