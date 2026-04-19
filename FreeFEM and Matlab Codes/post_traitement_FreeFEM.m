
clear; clc; close all;

filename = 'donnees_simulation_maillage_bon.txt';

if exist(filename, 'file')
    data = readmatrix(filename, 'Delimiter', ';', 'NumHeaderLines', 1);
    
    x = data(:, 1); % Colonne n
    y = data(:, 2); % Colonne Energie
    
    figure;
    plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    title('Évolution de l''Énergie en fonction de la densité de maillage');
    xlabel('Densite n');
    ylabel('Energie calculée (\int |\nabla w|^2)');
    grid on;
    
    fprintf('Nombre de points tracés : %d\n', length(x));
    fprintf('Energie min : %.4f\n', min(y));
    fprintf('Energie max : %.4f\n', max(y));
    
else
    error('Le fichier %s est introuvable. Vérifie le dossier courant.', filename);
end