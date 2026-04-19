clear; clc; close all;

% --- PARAMÈTRES ---
n_points = 10000;
eps_vals = logspace(-10, 0, n_points);
a_vals   = logspace(-10, 0, n_points); 
K = 100;    % Facteur de condition Z2 > K * Z1
n = 100;    % Paramètre nin (densité sur le trou)
N_max = 2000; % Limite pour éviter l'erreur de maillage bamg

% --- CALCULS SUR LA GRILLE ---
[E, A] = meshgrid(eps_vals, a_vals);

% Fonctions
f1 = @(e) -e.^2 .* log(e);
f2 = @(e,a) -e.^2 .* log(a);
Z1 = f1(E);
Z2 = -E.^2 .* log(A); % Calcul basé sur C0

% --- MASQUES DE SÉLECTION ---
% 1. Condition mathématique : Z2 > K * Z1
mask_math = Z2 > K*Z1;

% 2. Condition numérique : a >= 10^-5
mask_limit = A >= 1e-5;

% 3. Condition de maillage : n * epsilon / a < 2000
% Cette condition évite l'erreur "Assertion fail : (kkkk++ < 2000)"
mask_mesh = (n .* E ./ A) < N_max;

% Masque Global (Intersection de toutes les conditions)
mask_valid = mask_math & mask_limit & mask_mesh;

% --- EXTRACTION DES TRIPLETS ---
eps_final = E(mask_valid);
a_final   = A(mask_valid);

% Calcul de C0 = epsilon^2 * ln(a_epsilon)
C0_final = (eps_final.^2) .* log(a_final);

% Construction du tableau [Epsilon, a_epsilon, C0]
Resultats = [eps_final, a_final, C0_final];

% --- ÉCHANTILLONNAGE DES RÉSULTATS ---
pas = 1000; % On prend un point tous les 50 points trouvés
if size(Resultats, 1) > pas
    Resultats_echantillon = Resultats(1:pas:end, :);
else
    Resultats_echantillon = Resultats; % Pas assez de points pour échantillonner
end

% --- AFFICHAGE ---
fprintf('--- ANALYSE NUMÉRIQUE AVEC ÉCHANTILLONNAGE ---\n');
fprintf('Total points valides : %d\n', size(Resultats, 1));
fprintf('Nombre de points affichés (pas = %d) : %d\n\n', pas, size(Resultats_echantillon, 1));

if ~isempty(Resultats_echantillon)
    % Création de la table pour un affichage élégant
    T = array2table(Resultats_echantillon, 'VariableNames', {'Epsilon', 'a_epsilon', 'C0'});
    
    % Optionnel : on rajoute le calcul de nout pour vérifier la condition de maillage
    T.Nout_FreeFEM = ceil(n .* T.Epsilon ./ T.a_epsilon);
    
    disp(T);
else
    fprintf('AUCUN POINT TROUVÉ : Les contraintes sont trop restrictives.\n');
end

%% 
% --- REPRÉSENTATION 3D ---
% On réduit la densité pour le tracé 3D afin d'éviter de saturer la mémoire
n_plot = 200; 
eps_plot = logspace(-10, 0, n_plot);
a_plot   = logspace(-10, 0, n_plot);
[E_p, A_p] = meshgrid(eps_plot, a_plot);

Z1_p = -E_p.^2 .* log(E_p);
Z2_p = -E_p.^2 .* log(A_p);

figure('Color', 'w', 'Position', [100, 100, 1000, 600]);

% Surface f2 (C0)
surf(E_p, A_p, Z2_p, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'interp');
hold on;

% Surface f1 (Référence)
surf(E_p, A_p, Z1_p, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'FaceColor', [0.2 0.2 0.2]);

% Habillage
set(gca, 'XScale', 'log', 'YScale', 'log'); % Échelles logarithmiques
xlabel('\epsilon');
ylabel('a_\epsilon');
zlabel('Valeur des fonctions');
legend('f_2 = -\epsilon^2 ln(a_\epsilon)', 'f_1 = -\epsilon^2 ln(\epsilon)', 'Location', 'best');
view(135, 30); % Ajuste l'angle de vue
grid on;



%%

clear; clc; close all;

% --- DONNÉES 1 : Comparaison à la valeur limite (u_0) ---
h1 = [0.114154, 0.0412481, 0.0281514, 0.0202519, 0.0171872, 0.0140144, 0.0118296, 0.0102925, 0.00907024, 0.00819744];
tri1 = [1164, 10253, 27622, 56141, 92336, 138625, 187240, 250548, 329345, 401226];
err1 = [2.42849, 0.546328, 0.839405, 0.924243, 0.960011, 0.9784, 0.989089, 0.995846, 1.00039, 1.00358];

% --- DONNÉES 2 : Comparaison à la valeur théorique raffinée (u_epsilon) ---
err2 = [3.23504, 0.284808, 0.0058462, 0.0899836, 0.125456, 0.143692, 0.154294, 0.160994, 0.165497, 0.168668];

% --- CRÉATION DE LA FIGURE ---
figure('Color', 'w', 'Position', [100, 100, 900, 800]);

% --- LIGNE 1 : COMPARAISON À LA LIMITE ASYMPTOTIQUE ---

% En haut à gauche (Haut-Gauche)
subplot(2, 2, 1);
loglog(h1, err1, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
grid on; xlabel('Taille du maillage h'); ylabel('Erreur relative (%)');
title('Convergence (u_0) : Erreur vs h');

% En haut à droite (Haut-Droite)
subplot(2, 2, 2);
loglog(tri1, err1, '-sb', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; xlabel('Nombre de triangles'); ylabel('Erreur relative (%)');
title('Efficacité (u_0) : Erreur vs Triangles');

% --- LIGNE 2 : COMPARAISON À LA VALEUR THÉORIQUE RAFFINÉE ---

% En bas à gauche (Bas-Gauche)
subplot(2, 2, 3);
loglog(h1, err2, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.6 1]);
grid on; xlabel('Taille du maillage h'); ylabel('Erreur relative (%)');
title('Convergence (u_\epsilon) : Erreur vs h');

% En bas à droite (Bas-Droite)
subplot(2, 2, 4);
loglog(tri1, err2, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.4 0.4]);
grid on; xlabel('Nombre de triangles'); ylabel('Erreur relative (%)');
title('Efficacité (u_\epsilon) : Erreur vs Triangles');


