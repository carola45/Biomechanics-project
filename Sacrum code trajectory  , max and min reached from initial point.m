clear;
clc;
close;

% ====== PARAMETRI ======
fileName1 = 'Trial1.Anita.xlsx';  % Nome del primo file Excel
fileName2 = 'Trial2.Anita.xlsx';  % Nome del secondo file Excel
fileName3 = 'Trial3.Anita.xlsx';  % Nome del terzo file Excel
fileName4 = 'Trial4.Anita.xlsx';  % Nome del quarto file Excel
fileName5 = 'Trial5.Anita.xlsx';  % Nome del quinto file Excel
sheetName = 'Sheet1';  % Nome del foglio (se necessario)

% ====== CARICAMENTO DEI DATI ======
% Carica i dati dal primo file
data1 = readtable(fileName1, 'Sheet', sheetName, 'Range', 'A12');
X1 = data1{:, 15};  % Colonna 15 (X del primo file)
Y1 = data1{:, 16};  % Colonna 16 (Y del primo file)
Z1 = data1{:, 17};  % Colonna 17 (Z del primo file)

% Carica i dati dal secondo file
data2 = readtable(fileName2, 'Sheet', sheetName, 'Range', 'A12');
X2 = data2{:, 15};  % Colonna 15 (X del secondo file)
Y2 = data2{:, 16};  % Colonna 16 (Y del secondo file)
Z2 = data2{:, 17};  % Colonna 17 (Z del secondo file)

% Carica i dati dal terzo file
data3 = readtable(fileName3, 'Sheet', sheetName, 'Range', 'A12');
X3 = data3{:, 15};  % Colonna 15 (X del terzo file)
Y3 = data3{:, 16};  % Colonna 16 (Y del terzo file)
Z3 = data3{:, 17};  % Colonna 17 (Z del terzo file)

% Carica i dati dal quarto file
data4 = readtable(fileName4, 'Sheet', sheetName, 'Range', 'A12');
X4 = data4{:, 15};  % Colonna 15 (X del quarto file)
Y4 = data4{:, 16};  % Colonna 16 (Y del quarto file)
Z4 = data4{:, 17};  % Colonna 17 (Z del quarto file)

% Carica i dati dal quinto file
data5 = readtable(fileName5, 'Sheet', sheetName, 'Range', 'A12');
X5 = data5{:, 15};  % Colonna 15 (X del quinto file)
Y5 = data5{:, 16};  % Colonna 16 (Y del quinto file)
Z5 = data5{:, 17};  % Colonna 17 (Z del quinto file)

% ====== VISUALIZZAZIONE DELLE TRAIETTORIE NEL PIANO XY (SOTTOGRAPHI SEPARATI) ======
figure;

% Prima traiettoria (Trial 1) in rosso
subplot(5, 1, 1);  % 5 righe, 1 colonna, 1° grafico
plot(X1, Y1, '-o', 'LineWidth', 2, 'Color', 'r');
title('Traiettoria nel Piano XY - Trial 1');
xlabel('X');
ylabel('Y');
grid on;
xlim([min(X1)-1, max(X1)+1]);
ylim([min(Y1)-1, max(Y1)+1]);

% Seconda traiettoria (Trial 2) in verde
subplot(5, 1, 2);  % 5 righe, 1 colonna, 2° grafico
plot(X2, Y2, '-x', 'LineWidth', 2, 'Color', 'g');
title('Traiettoria nel Piano XY - Trial 2');
xlabel('X');
ylabel('Y');
grid on;
xlim([min(X2)-1, max(X2)+1]);
ylim([min(Y2)-1, max(Y2)+1]);

% Terza traiettoria (Trial 3) in blu
subplot(5, 1, 3);  % 5 righe, 1 colonna, 3° grafico
plot(X3, Y3, '-s', 'LineWidth', 2, 'Color', 'b');
title('Traiettoria nel Piano XY - Trial 3');
xlabel('X');
ylabel('Y');
grid on;
xlim([min(X3)-1, max(X3)+1]);
ylim([min(Y3)-1, max(Y3)+1]);

% Quarta traiettoria (Trial 4) in magenta
subplot(5, 1, 4);  % 5 righe, 1 colonna, 4° grafico
plot(X4, Y4, '-^', 'LineWidth', 2, 'Color', 'm');
title('Traiettoria nel Piano XY - Trial 4');
xlabel('X');
ylabel('Y');
grid on;
xlim([min(X4)-1, max(X4)+1]);
ylim([min(Y4)-1, max(Y4)+1]);

% Quinta traiettoria (Trial 5) in ciano
subplot(5, 1, 5);  % 5 righe, 1 colonna, 5° grafico
plot(X5, Y5, '-d', 'LineWidth', 2, 'Color', 'c');
title('Traiettoria nel Piano XY - Trial 5');
xlabel('X');
ylabel('Y');
grid on;
xlim([min(X5)-1, max(X5)+1]);
ylim([min(Y5)-1, max(Y5)+1]);

%% ====== VISUALIZZAZIONE DELLE TRAIETTORIE NEL PIANO ZY (SOTTOGRAPHI SEPARATI) ======
figure;

% Prima traiettoria (Trial 1) nel piano ZY in rosso
subplot(5, 1, 1);  % 5 righe, 1 colonna, 1° grafico
plot(Z1, Y1, '-o', 'LineWidth', 2, 'Color', 'r');
title('Traiettoria nel Piano ZY - Trial 1');
xlabel('Z');
ylabel('Y');
grid on;
ylim([min(Y1)-1, max(Y1)+1]);
xlim([min(Z1)-1, max(Z1)+1]);

% Seconda traiettoria (Trial 2) nel piano ZY in verde
subplot(5, 1, 2);  % 5 righe, 1 colonna, 2° grafico
plot(Z2, Y2, '-x', 'LineWidth', 2, 'Color', 'g');
title('Traiettoria nel Piano ZY - Trial 2');
xlabel('Z');
ylabel('Y');
grid on;
xlim([min(Z2)-1, max(Z2)+1]);
ylim([min(Y2)-1, max(Y2)+1]);

% Terza traiettoria (Trial 3) nel piano ZY in blu
subplot(5, 1, 3);  % 5 righe, 1 colonna, 3° grafico
plot(Z3, Y3, '-s', 'LineWidth', 2, 'Color', 'b');
title('Traiettoria nel Piano ZY - Trial 3');
xlabel('Z');
ylabel('Y');
grid on;
xlim([min(Z3)-1, max(Z3)+1]);
ylim([min(Y3)-1, max(Y3)+1]);

% Quarta traiettoria (Trial 4) nel piano ZY in magenta
subplot(5, 1, 4);  % 5 righe, 1 colonna, 4° grafico
plot(Z4, Y4, '-^', 'LineWidth', 2, 'Color', 'm');
title('Traiettoria nel Piano ZY - Trial 4');
xlabel('Z');
ylabel('Y');
grid on;
xlim([min(Z4)-1, max(Z4)+1]);
ylim([min(Y4)-1, max(Y4)+1]);

% Quinta traiettoria (Trial 5) nel piano ZY in ciano
subplot(5, 1, 5);  % 5 righe, 1 colonna, 5° grafico
plot(Z5, Y5, '-d', 'LineWidth', 2, 'Color', 'c');
title('Traiettoria nel Piano ZY - Trial 5');
xlabel('Z');
ylabel('Y');
grid on;
xlim([min(Z5)-1, max(Z5)+1]);
ylim([min(Y5)-1, max(Y5)+1]);

% ====== PARAMETRI ======
fileName1 = 'Trial1.Anita.xlsx';  % Nome del primo file Excel
fileName2 = 'Trial2.Anita.xlsx';  % Nome del secondo file Excel
fileName3 = 'Trial3.Anita.xlsx';  % Nome del terzo file Excel
fileName4 = 'Trial4.Anita.xlsx';  % Nome del quarto file Excel
fileName5 = 'Trial5.Anita.xlsx';  % Nome del quinto file Excel
sheetName = 'Sheet1';  % Nome del foglio (se necessario)

% ====== CARICAMENTO DEI DATI ======
% Carica i dati dal primo file
data1 = readtable(fileName1, 'Sheet', sheetName, 'Range', 'A12'); % Inizia dalla riga 12
X1 = data1{:, 15};  % Colonna 15 (X del primo file)
Y1 = data1{:, 16};  % Colonna 16 (Y del primo file)
Z1 = data1{:, 17};  % Colonna 17 (Z del primo file)

% Carica i dati dal secondo file
data2 = readtable(fileName2, 'Sheet', sheetName, 'Range', 'A12');
X2 = data2{:, 15};  % Colonna 15 (X del secondo file)
Y2 = data2{:, 16};  % Colonna 16 (Y del secondo file)
Z2 = data2{:, 17};  % Colonna 17 (Z del secondo file)

% Carica i dati dal terzo file
data3 = readtable(fileName3, 'Sheet', sheetName, 'Range', 'A12');
X3 = data3{:, 15};  % Colonna 15 (X del terzo file)
Y3 = data3{:, 16};  % Colonna 16 (Y del terzo file)
Z3 = data3{:, 17};  % Colonna 17 (Z del terzo file)

% Carica i dati dal quarto file
data4 = readtable(fileName4, 'Sheet', sheetName, 'Range', 'A12');
X4 = data4{:, 15};  % Colonna 15 (X del quarto file)
Y4 = data4{:, 16};  % Colonna 16 (Y del quarto file)
Z4 = data4{:, 17};  % Colonna 17 (Z del quarto file)

% Carica i dati dal quinto file
data5 = readtable(fileName5, 'Sheet', sheetName, 'Range', 'A12');
X5 = data5{:, 15};  % Colonna 15 (X del quinto file)
Y5 = data5{:, 16};  % Colonna 16 (Y del quinto file)
Z5 = data5{:, 17};  % Colonna 17 (Z del quinto file)

%% ====== CALCOLO DELLA VARIAZIONE MASSIMA E MINIMA RISPETTO AL VALORE INIZIALE (Y) ======
% Normalizza i valori di Y per ogni trial rispetto al primo valore (assunto come 0)
Y1_norm = Y1 - Y1(1);  % Normalizza Y1 rispetto al primo valore
Y2_norm = Y2 - Y2(1);  % Normalizza Y2 rispetto al primo valore
Y3_norm = Y3 - Y3(1);  % Normalizza Y3 rispetto al primo valore
Y4_norm = Y4 - Y4(1);  % Normalizza Y4 rispetto al primo valore
Y5_norm = Y5 - Y5(1);  % Normalizza Y5 rispetto al primo valore

% Calcola il massimo e il minimo di Y normalizzato per ogni trial
[maxY1, idxMaxY1] = max(Y1_norm);  % Massimo Y per il primo trial
[minY1, idxMinY1] = min(Y1_norm);  % Minimo Y per il primo trial

[maxY2, idxMaxY2] = max(Y2_norm);  % Massimo Y per il secondo trial
[minY2, idxMinY2] = min(Y2_norm);  % Minimo Y per il secondo trial

[maxY3, idxMaxY3] = max(Y3_norm);  % Massimo Y per il terzo trial
[minY3, idxMinY3] = min(Y3_norm);  % Minimo Y per il terzo trial

[maxY4, idxMaxY4] = max(Y4_norm);  % Massimo Y per il quarto trial
[minY4, idxMinY4] = min(Y4_norm);  % Minimo Y per il quarto trial

[maxY5, idxMaxY5] = max(Y5_norm);  % Massimo Y per il quinto trial
[minY5, idxMinY5] = min(Y5_norm);  % Minimo Y per il quinto trial

% ====== STAMPA DEI RISULTATI ======
disp('Trial 1:');
disp(['Max Y (normalizzato) = ', num2str(maxY1), ' alla riga ', num2str(idxMaxY1)]);
disp(['Min Y (normalizzato) = ', num2str(minY1), ' alla riga ', num2str(idxMinY1)]);

disp('Trial 2:');
disp(['Max Y (normalizzato) = ', num2str(maxY2), ' alla riga ', num2str(idxMaxY2)]);
disp(['Min Y (normalizzato) = ', num2str(minY2), ' alla riga ', num2str(idxMinY2)]);

disp('Trial 3:');
disp(['Max Y (normalizzato) = ', num2str(maxY3), ' alla riga ', num2str(idxMaxY3)]);
disp(['Min Y (normalizzato) = ', num2str(minY3), ' alla riga ', num2str(idxMinY3)]);

disp('Trial 4:');
disp(['Max Y (normalizzato) = ', num2str(maxY4), ' alla riga ', num2str(idxMaxY4)]);
disp(['Min Y (normalizzato) = ', num2str(minY4), ' alla riga ', num2str(idxMinY4)]);

disp('Trial 5:');
disp(['Max Y (normalizzato) = ', num2str(maxY5), ' alla riga ', num2str(idxMaxY5)]);
disp(['Min Y (normalizzato) = ', num2str(minY5), ' alla riga ', num2str(idxMinY5)]);

%% ====== CALCOLO DELLA VARIAZIONE MASSIMA E MINIMA RISPETTO AL VALORE INIZIALE ======
% Trova il massimo e minimo lungo l'asse Y (normalizzato) rispetto al valore iniziale
[maxY1, idxMaxY1] = max(Y1_norm);  % Massimo Y per il primo trial
[minY1, idxMinY1] = min(Y1_norm);  % Minimo Y per il primo trial

[maxY2, idxMaxY2] = max(Y2_norm);  % Massimo Y per il secondo trial
[minY2, idxMinY2] = min(Y2_norm);  % Minimo Y per il secondo trial

[maxY3, idxMaxY3] = max(Y3_norm);  % Massimo Y per il terzo trial
[minY3, idxMinY3] = min(Y3_norm);  % Minimo Y per il terzo trial

[maxY4, idxMaxY4] = max(Y4_norm);  % Massimo Y per il quarto trial
[minY4, idxMinY4] = min(Y4_norm);  % Minimo Y per il quarto trial

[maxY5, idxMaxY5] = max(Y5_norm);  % Massimo Y per il quinto trial
[minY5, idxMinY5] = min(Y5_norm);  % Minimo Y per il quinto trial

% Trova i valori corrispondenti di Z per i valori massimi e minimi di Y
maxZ1 = Z1(idxMaxY1);
minZ1 = Z1(idxMinY1);

maxZ2 = Z2(idxMaxY2);
minZ2 = Z2(idxMinY2);

maxZ3 = Z3(idxMaxY3);
minZ3 = Z3(idxMinY3);

maxZ4 = Z4(idxMaxY4);
minZ4 = Z4(idxMinY4);

maxZ5 = Z5(idxMaxY5);
minZ5 = Z5(idxMinY5);

% Calcola la variazione massima in Y rispetto alla posizione iniziale
varMaxY1 = maxY1;  % Variazione massima per il primo trial
varMinY1 = minY1;  % Variazione minima per il primo trial

varMaxY2 = maxY2;  % Variazione massima per il secondo trial
varMinY2 = minY2;  % Variazione minima per il secondo trial

varMaxY3 = maxY3;  % Variazione massima per il terzo trial
varMinY3 = minY3;  % Variazione minima per il terzo trial

varMaxY4 = maxY4;  % Variazione massima per il quarto trial
varMinY4 = minY4;  % Variazione minima per il quarto trial

varMaxY5 = maxY5;  % Variazione massima per il quinto trial
varMinY5 = minY5;  % Variazione minima per il quinto trial

% Stampa i risultati
disp('Trial 1:');
disp(['Max Y (normalized) = ', num2str(varMaxY1), ' at Z = ', num2str(maxZ1)]);
disp(['Min Y (normalized) = ', num2str(varMinY1), ' at Z = ', num2str(minZ1)]);

disp('Trial 2:');
disp(['Max Y (normalized) = ', num2str(varMaxY2), ' at Z = ', num2str(maxZ2)]);
disp(['Min Y (normalized) = ', num2str(varMinY2), ' at Z = ', num2str(minZ2)]);

disp('Trial 3:');
disp(['Max Y (normalized) = ', num2str(varMaxY3), ' at Z = ', num2str(maxZ3)]);
disp(['Min Y (normalized) = ', num2str(varMinY3), ' at Z = ', num2str(minZ3)]);

disp('Trial 4:');
disp(['Max Y (normalized) = ', num2str(varMaxY4), ' at Z = ', num2str(maxZ4)]);
disp(['Min Y (normalized) = ', num2str(varMinY4), ' at Z = ', num2str(minZ4)]);

disp('Trial 5:');
disp(['Max Y (normalized) = ', num2str(varMaxY5), ' at Z = ', num2str(maxZ5)]);
disp(['Min Y (normalized) = ', num2str(varMinY5), ' at Z = ', num2str(minZ5)]);
