close all
clear all 
clc
set(groot, ...
    'defaultAxesFontSize', 28, ...
    'defaultLegendFontSize', 12);   % keep legend smaller

output_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon\V14/my_post_proc_figures\Jint\';
mat_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon\V14/my_post_proc_figures\Jint\postproc Jint fit results.mat';
postStructmat_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon\V14/my_post_proc_figures\postProcessingStruct_Intrinsic.mat';  % <-- Your custom folder path

Jcap1 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
Jcap2 = [0.15, 0.135, 0.12, 0.105, 0.09, 0.076, 0.061, 0.046, 0.031, 0.016, 0.005, 0.001];

filter_AB = false;
min_PSI_points = 10;

run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB);
plot_fit_results_from_TXT(mat_path, output_path);

function run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB)

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    if isempty(Jcap1) || isempty(Jcap2)
        error('Jcap1 and Jcap2 must not be empty');
    end

    processedPairs = containers.Map();  % keys like '0.001_0.150'

    for j = 1:length(Jcap1)
        for k = 1:length(Jcap2)
            N1 = Jcap1(j);
            N2 = Jcap2(k);

            if abs(N1 - N2) < 1e-10
                fprintf('Skipping identical pair: [%.4f, %.4f]\n', N1, N2);
                continue;
            end

            % Create a unique key for unordered pair
            key = sprintf('%.10f_%.10f', min(N1, N2), max(N1, N2));
            if isKey(processedPairs, key)
                fprintf('Skipping mirrored duplicate: [%.4f, %.4f]\n', N1, N2);
                continue;
            end

            % Mark as processed
            processedPairs(key) = true;

            fprintf('Processing pair [%.4f, %.4f]\n', N1, N2);
            try
                postproc_fit_results(N1, N2, postStructmat_path, output_path, filter_AB, min_PSI_points);
            catch ME
                warning('Failed processing pair [%.4f, %.4f]: %s', N1, N2, ME.message);
            end
        end
    end
end

function postproc_fit_results(Jcap1, Jcap2, postStructmat_path, output_path, filter_AB, min_PSI_points)
% POSTPROC_FIT_RESULTS
% 1) Opens the .mat struct and collects data keyed by mismatch pair [a,b].
% 2) For each mismatch pair [a,b] with a ~= b, the function also looks up:
%       [a,a] and [b,b] (if they exist in the struct),
%    combines those data, performs the fit:
%
%       Z = ((psi*(1-psi))^2)*(phi^a)*b
%
%    and plots all points in one figure (with different markers for each pair).
% 3) Writes results to "postproc Jint fit results.txt" and saves each figure in the same folder.

    %% 1) Load struct and find relevant data
    S = load(postStructmat_path);
    fn = fieldnames(S);
    if numel(fn) ~= 1 || ~isstruct(S.(fn{1}))
        error('File %s must contain ONE struct variable.', postStructmat_path);
    end
    S = S.(fn{1});

    tol = 1e-6;
    aVal = Jcap1;
    bVal = Jcap2;
    if abs(aVal - bVal) < tol
        error('Jcap1 and Jcap2 must be different.');
    end

    % Preallocate
    phiAB = []; psiAB = []; Z_AB = [];
    dataAA = struct('phi',[],'psi',[],'Z',[]);
    dataBB = dataAA;

    % Extract AB, AA, BB data
    for k = 1:numel(S)
        jcaps = unique(S(k).Jcap(:));
        jcaps = sort(jcaps);
        if numel(jcaps)==2 && ismembertol(aVal,jcaps,tol) && ismembertol(bVal,jcaps,tol)
            phiAB = [phiAB; S(k).phi(:)];
            psiAB = [psiAB; S(k).psi(:)];
            Z_AB  = [Z_AB ; S(k).Jintrin(:)];
            continue
        end
        if numel(jcaps)==1 && abs(jcaps - aVal) < tol
            dataAA.phi = [dataAA.phi; S(k).phi(:)];
            dataAA.Z   = [dataAA.Z  ; S(k).Jintrin(:)];
            continue
        end
        if numel(jcaps)==1 && abs(jcaps - bVal) < tol
            dataBB.phi = [dataBB.phi; S(k).phi(:)];
            dataBB.Z   = [dataBB.Z  ; S(k).Jintrin(:)];
        end
    end

    if isempty(phiAB)
        error('Pair [%.4f, %.4f] not found in %s', aVal, bVal, postStructmat_path);
    end

    phiAB_raw = phiAB; psiAB_raw = psiAB; Z_AB_raw = Z_AB;

    if filter_AB
        tolPsi = 1e-6;
        minRSquare = 0.9;
        minPointsPerPsi = min_PSI_points;

        [uniquePsi, ~, idPsi] = uniquetol(psiAB, tolPsi);
        psiCounts = accumarray(idPsi, 1);
        validBands = psiCounts(idPsi) >= minPointsPerPsi;

        phiAB = phiAB(validBands);
        psiAB = psiAB(validBands);
        Z_AB  = Z_AB(validBands);

        [phiAB_filt, psiAB_filt, Z_AB_filt] = ...
            filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, tolPsi, minRSquare, minPointsPerPsi);

        if isempty(phiAB_filt)
            warning('AB data is empty after filtering.');
            return;
        end
    else
        phiAB_filt = phiAB;
        psiAB_filt = psiAB;
        Z_AB_filt  = Z_AB;
    end

    if aVal > bVal
        dataAA.psi = zeros(size(dataAA.phi));
        dataBB.psi = ones(size(dataBB.phi));
    else
        dataAA.psi = ones(size(dataAA.phi));
        dataBB.psi = zeros(size(dataBB.phi));
    end

    psiAll = [psiAB_filt; dataAA.psi; dataBB.psi];
    phiAll = [phiAB_filt; dataAA.phi; dataBB.phi];
    ZAll   = [Z_AB_filt ; dataAA.Z   ; dataBB.Z];

    if isempty(psiAll), return; end

    [fitresult, gof, psiAB_in, phiAB_in, ZAB_in] = ...
        fitWithMultipleStartPoints(aVal, bVal, psiAll, phiAll, ZAll);

    aFitted = fitresult.a;
    bFitted = fitresult.b;
    r2      = gof.rsquare;

    ci = confint(fitresult, 0.95);
    aErrPlus  = ci(2,1) - aFitted;
    aErrMinus = aFitted - ci(1,1);
    bErrPlus  = ci(2,2) - bFitted;
    bErrMinus = bFitted - ci(1,2);
    a_covers = (aVal >= ci(1,1)) && (aVal <= ci(2,1));
    b_covers = (bVal >= ci(1,2)) && (bVal <= ci(2,2));

%% ------------------------------------------------------------------------
%  build / append the per-fit table row and save one MAT file
%% ------------------------------------------------------------------------
[outFolder,~,~] = fileparts(output_path);
if ~exist(outFolder,'dir'), mkdir(outFolder); end

matFile = fullfile(outFolder,'postproc Jint fit results.mat');

% ---------- persistent guard: wipe file only once per MATLAB session -----
persistent clearedOnce
if isempty(clearedOnce)
    if isfile(matFile), delete(matFile); end
    clearedOnce = true;
end

% ---------- decide if we have a c coefficient ---------------------------
haveC = exist('cFit','var') && isscalar(cFit);

if haveC
    colNames = {'Jcap1','Jcap2', ...
                'a','b','c', ...
                'daPlus','daMinus','dbPlus','dbMinus', ...
                'Rsq'};
else
    colNames = {'Jcap1','Jcap2', ...
                'a','b', ...
                'daPlus','daMinus','dbPlus','dbMinus', ...
                'Rsq'};
end

% ---------- load existing table or start a new one ----------------------
if isfile(matFile)
    S = load(matFile);
    if isfield(S,'fitTable')
        fitTable = S.fitTable;
    else
        fitTable = table('Size',[0 numel(colNames)], ...
                         'VariableTypes', repmat("double",1,numel(colNames)), ...
                         'VariableNames', colNames);
    end
else
    fitTable = table('Size',[0 numel(colNames)], ...
                     'VariableTypes', repmat("double",1,numel(colNames)), ...
                     'VariableNames', colNames);
end

% ---------- ensure all scalars ------------------------------------------
aErrPlus   = aErrPlus(1);
aErrMinus  = aErrMinus(1);
bErrPlus   = bErrPlus(1);
bErrMinus  = bErrMinus(1);
if ~haveC, cFit = []; end          % so we can build the row generically

% ---------- build one-row table -----------------------------------------
if haveC
    newRow = table( aVal     , bVal , ...
                    aFitted  , bFitted , cFit , ...
                    aErrPlus , aErrMinus , ...
                    bErrPlus , bErrMinus , ...
                    r2 , ...
                    'VariableNames', colNames );
else
    newRow = table( aVal     , bVal , ...
                    aFitted  , bFitted , ...
                    aErrPlus , aErrMinus , ...
                    bErrPlus , bErrMinus , ...
                    r2 , ...
                    'VariableNames', colNames );
end

% ---------- append and save ---------------------------------------------
fitTable = [fitTable ; newRow];
equationStr = "Fit Equation: Jint = \phi[ \psi·J_a + (1-\psi)·J_b ]";
save(matFile,'fitTable','equationStr','-v7');   % use -v7.3 if >2 GB

fprintf('Fit results appended to MAT-file:\n   %s\n', matFile);



%% Plot: X=psi, Y=phi, Z=Jint
hFig = figure('Visible','on'); hold on; grid on; view(3);

% --- Raw AB
scatter3(psiAB_raw, phiAB_raw, Z_AB_raw, ...
    36,'d','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','none', ...
    'DisplayName',sprintf('[%.4f, %.4f]  AB - raw',aVal,bVal));

% --- Filtered AB
scatter3(psiAB_filt, phiAB_filt, Z_AB_filt, ...
    40,'d','MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor','none', ...
    'DisplayName',sprintf('[%.4f, %.4f]  AB - filtered',aVal,bVal));

% --- Inliers actually used in the fit (from fitWithMultipleStartPoints)
scatter3(psiAB_in, phiAB_in, ZAB_in, ...
    40,'d','MarkerEdgeColor','k','MarkerFaceColor','k', ...
    'DisplayName',sprintf('[%.4f, %.4f]  used',aVal,bVal));

% --- AA and BB points (psi=1 or 0 according to your convention)
if ~isempty(dataAA.phi)
    scatter3(dataAA.psi, dataAA.phi, dataAA.Z, ...
        60,'^','MarkerEdgeColor','r','MarkerFaceColor','r', ...
        'DisplayName',sprintf('[%.4f, %.4f]  AA',aVal,aVal));
end
if ~isempty(dataBB.phi)
    scatter3(dataBB.psi, dataBB.phi, dataBB.Z, ...
        60,'s','MarkerEdgeColor','b','MarkerFaceColor','b', ...
        'DisplayName',sprintf('[%.4f, %.4f]  BB',bVal,bVal));
end

% --- fitted surface: Z = phi*(psi*a + (1-psi)*b)
psiMin = 0; psiMax = 1;
phiMin = max(0, min(phiAll));
phiMax = max(phiAll);

[psi_grid, phi_grid] = meshgrid(linspace(psiMin, psiMax, 60), ...
                                linspace(phiMin, phiMax, 60));

Z_grid = phi_grid .* (psi_grid .* aFitted + (1 - psi_grid) .* bFitted);

surf(psi_grid, phi_grid, Z_grid, ...
    'EdgeColor','none','FaceAlpha',0.35, ...
    'DisplayName','Fitted surface');

xlabel('\psi','Interpreter','tex','FontSize',40,'Color','red','FontWeight','bold');
ylabel('\phi','Interpreter','tex','FontSize',40,'Color','red','FontWeight','bold');
zlabel('$J_{\mathrm{int}}\;\left[\frac{1}{nm}\right]$', ...
       'Interpreter','latex','FontSize',40,'Color','red','FontWeight','bold');

title(sprintf('J_{int} fit   [J_a , J_b] = [%.4f , %.4f]', aVal, bVal), ...
      'FontSize',40,'Color','red','FontWeight','bold');

eqStr = 'Fit surface:  J_{int}(\psi,\phi) = \phi[\psi J_{a} + (1-\psi)J_{b}]';
% fitStr = sprintf('J_{a} = %.4g,   J_{b} = %.4g,   R^2 = %.3f', aFitted, bFitted, r2);

subtitle({eqStr}, 'Interpreter','tex', ...
    'FontSize',26,'Color','red','FontWeight','bold');

legend('Location','best');
legend('Location','best','FontSize',12);

hold off;

savefig(hFig, fullfile(outFolder, sprintf('fit_Jint_vs_psi_phi_%.4f_%.4f.fig', aVal, bVal)));
% saveas(hFig, fullfile(outFolder, sprintf('fit_Jint_vs_psi_phi_%.4f_%.4f.png', aVal, bVal)));
close(hFig);


    fprintf('Results have been written to: %s\n', outFolder);
end


function plot_fit_results_from_TXT(txtFilePath, output_path)
% PLOT_FIT_RESULTS_FROM_TXT Reads the "postproc fit results.txt" file and creates
% four 2D figures:
%   1) a values vs. mismatch delta (delta = highVal - lowVal).
%   2) b values vs. mismatch delta.
%   3) R^2 values vs. mismatch delta.
%   4) Ratios of curvatures: 
%       min(Ja,Jb)/min(a,b)  and  max(Ja,Jb)/max(a,b)  vs. mismatch delta.
%
% The function expects text blocks like:
%
%     Fit Equation: Fmin = ((psi*(1-psi))^2)*(phi^a)*b
%
%     Fitting pairs [0.0010, 0.0050] + [0.0010, 0.0010] + [0.0050, 0.0050]
%        a = 1.00001
%        b = 3.90752e-05
%        R^2 = -1.7037
%
% The first bracketed pair is the "actual" (Ja,Jb),
% the second is "low homogeneous," and the third is "high homogeneous".
%
% The mismatch delta is computed as: delta = (high homogeneous) - (low homogeneous).
%
% Usage:
%   plot_fit_results_from_TXT('C:\path\to\postproc fit results.txt');

    % Extract folder path from txtFilePath
    [folderPath, ~, ~] = fileparts(txtFilePath);
    outFileName = fullfile(folderPath, 'postproc fit results.txt');

% === NEW: read the table directly ====================================
S = load(txtFilePath);          % txtFilePath now points to the MAT file
T = S.fitTable;

deltaVals      = abs(T.Jcap2 - T.Jcap1);
aValues        = T.a;
bValues        = T.b;
r2Values       = T.Rsq;
jaVals         = T.Jcap1;
jbVals         = T.Jcap2;

lowVals        = min(T.Jcap1, T.Jcap2);
highVals       = max(T.Jcap1, T.Jcap2);
mismatchLabels = arrayfun(@(lo,hi)sprintf('[%.4g, %.4g]',lo,hi), ...
                          lowVals, highVals, 'UniformOutput',false);
% =====================================================================


    %% ================= SORT ALL DATA BY DELTA =================
    [sortedDelta, sortIdx] = sort(deltaVals);

    sortedA                   = aValues(sortIdx);
    sortedB                   = bValues(sortIdx);
    sortedR2                  = r2Values(sortIdx);
    sortedMismatchLabels      = mismatchLabels(sortIdx);
    sortedJaVals              = jaVals(sortIdx);
    sortedJbVals              = jbVals(sortIdx);

    %% =============== Create Figures ===============

    % --- Figure 1:  a vs. mismatch delta ------------------------------------
    hFig1 = figure('Visible','on');
    plot(sortedDelta, sortedA, 'rs-', 'LineWidth',2, 'MarkerSize',8);
    xlabel('Mismatch \Delta (High - Low)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    ylabel('a Value',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    title('a Values vs. Mismatch \Delta',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    subtitle('Intrinsic curvature model');
    grid on;
    
    %% >>> add point-labels <<< ----------------------------------------------
    for i = 1:numel(sortedDelta)
        text(sortedDelta(i), sortedA(i), sortedMismatchLabels{i}, ...
             'VerticalAlignment','bottom', 'HorizontalAlignment','right', ...
             'FontSize',8, 'Color','k');
    end
    % ------------------------------------------------------------------------
    
    savefig(hFig1, fullfile(folderPath, 'aValues_vs_Delta.fig'));
    saveas(hFig1, fullfile(folderPath, 'aValues_vs_Delta.png'));
    close(hFig1);

    % --- Figure 2: b vs. mismatch delta ---
    hFig2 = figure('Visible','on');
    plot(sortedDelta, sortedB, 'bo-', 'LineWidth',2, 'MarkerSize',8);
    xlabel('Mismatch \Delta (High - Low)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    ylabel('b Value',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    title('b Values vs. Mismatch \Delta',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    subtitle('Intrinsic curvature model');
    grid on;
    % Annotate each point with the mismatch label
    for i = 1:length(sortedDelta)
        text(sortedDelta(i), sortedB(i), sortedMismatchLabels{i}, ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','right', ...
            'FontSize', 8, 'Color', 'k');
    end
    savefig(hFig2, fullfile(folderPath, 'bValues_vs_Delta.fig'));
    saveas(hFig2, fullfile(folderPath, 'bValues_vs_Delta.png'));
    close(hFig2);

    % --- Figure 3: R^2 vs. mismatch delta ---
    hFig3 = figure('Visible','on');
    plot(sortedDelta, sortedR2, 'cd-', 'LineWidth',2, 'MarkerSize',8);
    xlabel('Mismatch \Delta (High - Low)',   'FontSize',40, 'Color','red', 'FontWeight','bold');
    ylabel('R^2',   'FontSize',40, 'Color','red', 'FontWeight','bold');
    title('R^2 Values vs. Mismatch \Delta',   'FontSize',40, 'Color','red', 'FontWeight','bold');
    subtitle('Intrinsic curvature model');
    grid on;

        %% >>> add point-labels <<< ----------------------------------------------
    for i = 1:numel(sortedDelta)
        text(sortedDelta(i), sortedR2(i), sortedMismatchLabels{i}, ...
             'VerticalAlignment','bottom', 'HorizontalAlignment','right', ...
             'FontSize',8, 'Color','k');
    end
    % ------------------------------------------------------------------------

    savefig(hFig3, fullfile(folderPath, 'R2Values_vs_Delta.fig'));
    saveas(hFig3, fullfile(folderPath, 'R2Values_vs_Delta.png'));
    close(hFig3);
    
   %% ---------- Figure 4a & 4b :  min-ratio  and  max-ratio  vs Δ ----------

denMin = min(sortedA , sortedB);
denMax = max(sortedA , sortedB);

% -- 1) compute the two ratios for every data point -----------------------
ratioMinVal = min(sortedJaVals , sortedJbVals) ./ max(denMin , eps);
ratioMaxVal = max(sortedJaVals , sortedJbVals) ./ max(denMax , eps);

% -- 2) optional cutoff to hide absurd points ----------------------------
cutoff      = 1e12;          % keep |ratio| < cutoff
validIdx = isfinite(ratioMinVal) & isfinite(ratioMaxVal);

ratioMinVal = ratioMinVal(validIdx);
ratioMaxVal = ratioMaxVal(validIdx);
deltaOK     = sortedDelta(validIdx);
jaOK        = sortedJaVals(validIdx);
jbOK        = sortedJbVals(validIdx);

% ---------- (a)  min-ratio figure ---------------------------------------
hFig4a = figure('Visible','on');   hold on;  grid on;
plot(deltaOK, ratioMinVal, 'r-o', ...
     'LineWidth',2,'MarkerSize',8, ...
     'DisplayName','min(J_a,J_b)/min(a,b)');

xlabel('Mismatch \Delta  (High – Low)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
ylabel('min(J_a , J_b)  /  min(a , b)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
title('Minimum-curvature ratio vs. mismatch \Delta',   'FontSize',30, 'Color','red', 'FontWeight','bold');
subtitle('Intrinsic curvature model');

% annotate each point with [J_a , J_b]
for i = 1:numel(deltaOK)
    txt = sprintf('[%.4g , %.4g]', jaOK(i), jbOK(i));
    text(deltaOK(i), ratioMinVal(i), txt, ...
         'VerticalAlignment','top', 'HorizontalAlignment','right', ...
         'FontSize',8, 'Color','r');
end
savefig(hFig4a, fullfile(folderPath,'ratioMin_vs_Delta.fig'));
saveas(hFig4a, fullfile(folderPath, 'ratioMin_vs_Delta.png'));
close(hFig4a);

% ---------- (b)  max-ratio figure ---------------------------------------
hFig4b = figure('Visible','on');   hold on;  grid on;
plot(deltaOK, ratioMaxVal, 'b-s', ...
     'LineWidth',2,'MarkerSize',8, ...
     'DisplayName','max(J_a,J_b)/max(a,b)');

xlabel('Mismatch \Delta  (High – Low)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
ylabel('max(J_a , J_b)  /  max(a , b)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
title('Maximum-curvature ratio vs. mismatch \Delta',   'FontSize',30, 'Color','red', 'FontWeight','bold');
subtitle('Intrinsic curvature model');

for i = 1:numel(deltaOK)
    txt = sprintf('[%.4g , %.4g]', jaOK(i), jbOK(i));
    text(deltaOK(i), ratioMaxVal(i), txt, ...
         'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
         'FontSize',8, 'Color','b');
end
savefig(hFig4b, fullfile(folderPath,'ratioMax_vs_Delta.fig'));
saveas(hFig4b, fullfile(folderPath, 'ratioMax_vs_Delta.png'));
close(hFig4b);


[sigma_a_total, sigma_b_total, da_all, db_all, ...
 Jcap_min, Jcap_max, ab_min, ab_max] = compute_overall_uncertainty(txtFilePath, output_path);

% --- Determine which value came from a or b ---------------------
isAmin = ab_min < ab_max;   % if ab_min is a → error is da_all
isAmax = ab_max > ab_min;   % if ab_max is b → error is db_all

fprintf("Size check → isAmin: %d, da_all: %d, db_all: %d\n", ...
    numel(isAmin), numel(da_all), numel(db_all));

dy_min = db_all;                    % default to db
dy_min(isAmin) = da_all(isAmin);   % assign da where appropriate

dy_max = da_all;                    % default to da
dy_max(isAmax) = db_all(isAmax);   % assign db where appropriate

% =======================
% Weighted linear fit #1 (min) with no intercept
% =======================
w_min = 1 ./ dy_min.^2;
w_min(w_min > 1e8) = 1e8;
w_min = w_min / max(w_min);

ftype = fittype('m * x', 'independent', 'x', 'coefficients', 'm');
mdl_min = fit(Jcap_min, ab_min, ftype, ...
    'Weights', w_min, 'StartPoint', 1, 'Lower', 0, 'Upper', 2);

yfit_min = mdl_min(Jcap_min);
ci_min = confint(mdl_min, 0.6827);  % 1σ
slope_se_min = 0.5 * diff(ci_min); % Only one parameter now

% Re-plot: ab_min vs Jcap_min with fit
fig1 = figure('Visible','on'); hold on; grid on;
errorbar(Jcap_min, ab_min, dy_min, ...
         'ko', 'LineWidth', 1.5, 'CapSize', 4, ...
         'MarkerFaceColor','k', 'DisplayName','min(a,b)');
plot(Jcap_min, yfit_min, 'r-', 'LineWidth', 2, 'DisplayName','Weighted fit');

SSE = sum((ab_min - yfit_min).^2);
SST = sum((ab_min - mean(ab_min)).^2);
R2_min  = 1 - SSE/SST;


xlabel('Actual J_{prot1}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
ylabel('Fitted J_{prot1}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
title('Actual J_{prot1}  vs  Fitted J_{prot1}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
eqText = sprintf('Fit:  fitted J_{prot1} = (%.4g \\pm %.1e)·Actual J_{prot1})    (R^2 = %.3f)', ...
                 mdl_min.m, slope_se_min, R2_min);
subtitle(eqText, 'Interpreter','tex', ...
         'FontSize',22, 'Color','red', 'FontWeight','bold');
legend('Location','best');
legend('Location','best','FontSize',12);

txt1 = sprintf('Y = [%.4f ± %.1e]·X', mdl_min.m, slope_se_min);
text(0.05, 0.92, txt1,'FontSize',40, 'Units','normalized', 'FontWeight','bold');

savefig(fig1, fullfile(folderPath, 'ab_min_vs_Jcap_min_with_fit.fig'));
saveas(fig1, fullfile(folderPath, 'ab_min_vs_Jcap_min_with_fit.png'));
close(fig1);

% =======================
% Weighted linear fit #2 (max) with no intercept
% =======================
w_max = 1 ./ dy_max.^2;
w_max(w_max > 1e8) = 1e8;
w_max = w_max / max(w_max);

% --- Print weights (optional debug) ---
fprintf('\n%-12s %-12s %-12s %-12s\n', 'Jcap_max', 'ab_max', 'σ_y', 'weight');
fprintf('-------------------------------------------------------------\n');
for i = 1:numel(Jcap_max)
    fprintf('%-12.5f %-12.5f %-12.3e %-12.3e\n', ...
            Jcap_max(i), ab_max(i), dy_max(i), w_max(i));
end
fprintf('\nWeight range: min = %.3e  |  max = %.3e\n', min(w_max), max(w_max));
fprintf('σ_y range:    min = %.3e  |  max = %.3e\n\n', min(dy_max), max(dy_max));

% --- Fit ---
mdl_max = fit(Jcap_max, ab_max, ftype, ...
    'Weights', w_max, 'StartPoint', 1, 'Lower', 0, 'Upper', 2);
yfit_max = mdl_max(Jcap_max);
ci_max = confint(mdl_max, 0.6827);
slope_se_max = 0.5 * diff(ci_max);

fig2 = figure('Visible','on'); hold on; grid on;
errorbar(Jcap_max, ab_max, dy_max, ...
         'bo', 'LineWidth', 1.5, 'CapSize', 4, ...
         'MarkerFaceColor','b', 'DisplayName','max(a,b)');
plot(Jcap_max, yfit_max, 'r-', 'LineWidth', 2, 'DisplayName','Weighted fit');

SSE = sum((ab_max - yfit_max).^2);
SST = sum((ab_max - mean(ab_max)).^2);
R2_max  = 1 - SSE/SST;

eq_max = sprintf('\\itmax\\rm(a,b) = (%.4g \\pm %.1e)\\cdot \\itmax\\rm(J_{cap1}, J_{cap2})', ...
                 mdl_max.m, slope_se_max);
r2_max = sprintf('R^2 = %.3f', R2_max);

xlabel('Actual J_{prot2}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
ylabel('Fitted J_{prot2}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
title('Actual J_{prot2}  vs  Fitted J_{prot2}',   'FontSize',40, 'Color','red', 'FontWeight','bold');
eqText = sprintf('Fit:  fitted J_{prot2} = (%.4g \\pm %.1e)·Actual J_{prot2})    (R^2 = %.3f)', ...
                 mdl_max.m, slope_se_max, R2_max);
subtitle(eqText, 'Interpreter','tex', ...
         'FontSize',22, 'Color','red', 'FontWeight','bold');
legend('Location','best');
legend('Location','best','FontSize',12);

txt2 = sprintf('Y = [%.4f ± %.1e]·X', mdl_max.m, slope_se_max);
text(0.05, 0.92, txt2, 'FontSize',40, 'Units','normalized', 'FontWeight','bold');

savefig(fig2, fullfile(folderPath, 'ab_max_vs_Jcap_max_with_fit.fig'));
saveas(fig2, fullfile(folderPath, 'ab_max_vs_Jcap_max_with_fit.png'));
close(fig2);

    fprintf('Results have been written to: %s\n', outFileName);
end



function [bestFit, bestGOF, psiClean, phiClean, ZClean] = ...
    fitWithMultipleStartPoints(~, ~, psiAll, phiAll, ZAll)
% Fit model:
%   Jint(psi,phi) = phi * ( psi*a + (1-psi)*b )
% where a,b are fitted coefficients.

    % columnize
    psiAll = psiAll(:);
    phiAll = phiAll(:);
    ZAll   = ZAll(:);

    % basic validity
    valid = isfinite(psiAll) & isfinite(phiAll) & isfinite(ZAll) & (phiAll > 0);
    psiAll = psiAll(valid);
    phiAll = phiAll(valid);
    ZAll   = ZAll(valid);

    psiClean = psiAll;
    phiClean = phiAll;
    ZClean   = ZAll;

    if numel(ZAll) < 10
        error('Not enough valid points for fitting (need >= 10).');
    end

    % === initial guesses ===
    % From:
    %   Z/phi = psi*a + (1-psi)*b = b + psi*(a-b)
    % So endpoints give quick guesses:
    y = ZAll ./ phiAll;

    % prefer bands near psi ~ 1 and psi ~ 0 if they exist
    near1 = psiAll > 0.9;
    near0 = psiAll < 0.1;

    if any(near1)
        a_guess = median(y(near1), 'omitnan');
    else
        a_guess = median(y, 'omitnan');
    end
    if any(near0)
        b_guess = median(y(near0), 'omitnan');
    else
        b_guess = median(y, 'omitnan');
    end

    if ~isfinite(a_guess), a_guess = 0.05; end
    if ~isfinite(b_guess), b_guess = 0.05; end

    startPoints = [
        a_guess,         b_guess
        0.5*a_guess,     0.5*b_guess
        2.0*a_guess,     2.0*b_guess
        1.2*a_guess,     0.8*b_guess
        0.8*a_guess,     1.2*b_guess
    ];

    % === fit type ===
    ft = fittype('phi*(psi*a + (1-psi)*b)', ...
        'independent', {'psi','phi'}, ...
        'coefficients', {'a','b'});

    opts = fitoptions(ft);
    opts.Robust      = 'Bisquare';
    opts.Lower       = [0, 0];
    opts.Upper       = [Inf, Inf];
    opts.MaxFunEvals = 1e4;
    opts.MaxIter     = 1e4;

    bestR2  = -Inf;
    bestFit = [];
    bestGOF = [];

    for i = 1:size(startPoints,1)
        opts.StartPoint = startPoints(i,:);
        try
            [fitTmp, gofTmp] = fit([psiAll, phiAll], ZAll, ft, opts);
            if gofTmp.rsquare > bestR2
                bestR2  = gofTmp.rsquare;
                bestFit = fitTmp;
                bestGOF = gofTmp;
            end
        catch
            continue
        end
    end

    if isempty(bestFit)
        error('All fits failed.');
    end
end



function [phiAB_filt, psiAB_filt, Z_AB_filt] = ...
    filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, ...
                                  tolPsi, minRSquare, minPointsPerPsi)

% FILTER_AB_BY_PSI_FIT_QUALITY (Linear fit version, simplified for postproc_fit_results)
% -------------------------------------------------------------------------
% Filters AB data by evaluating the fit quality (R²) of Z vs φ at each ψ slice.
% Discards entire ψ bands where linear fit R² is below minRSquare.
%
% Inputs:
%   phiAB, psiAB, Z_AB : AB data arrays
%   tolPsi        : tolerance for grouping ψ (default: 1e-6)
%   minRSquare    : minimum R² to keep a ψ slice (default: 0.95)
%   minPointsPerPsi : minimum points per ψ slice to attempt fitting (default: 5)
%
% Outputs:
%   Filtered phiAB, psiAB, Z_AB
% -------------------------------------------------------------------------

if nargin < 4, tolPsi = 1e-6; end
if nargin < 5, minRSquare = 0.95; end
if nargin < 6, minPointsPerPsi = 5; end

[uniquePsi,~,idPsi] = uniquetol(psiAB, tolPsi);

keepMask = false(size(psiAB));   % initialize as all false

for k = 1:numel(uniquePsi)
    idxThisPsi = idPsi == k;
    
    if sum(idxThisPsi) < minPointsPerPsi
        continue
    end
    
    % --- Linear fit ---
    phi_k = phiAB(idxThisPsi);
    Z_k   = Z_AB(idxThisPsi);
    
    p = polyfit(phi_k, Z_k, 1);      % linear fit
    Z_fit = polyval(p, phi_k);
    
    % R² calculation
    SS_res = sum((Z_k - Z_fit).^2);
    SS_tot = sum((Z_k - mean(Z_k)).^2);
    R2 = 1 - SS_res / SS_tot;
    
    % Keep only if R² good
    if R2 >= minRSquare
        keepMask(idxThisPsi) = true;
    end
end

% Apply mask to AB
phiAB_filt = phiAB(keepMask);
psiAB_filt = psiAB(keepMask);
Z_AB_filt  = Z_AB(keepMask);

end

function [sigma_a_total, sigma_b_total, da_all, db_all, ...
          Jcap_min, Jcap_max, ab_min, ab_max] = compute_overall_uncertainty(matFilePath, output_path)
% -------------------------------------------------------------------------
% COMPUTE_OVERALL_UNCERTAINTY
% Parses the result TXT file and returns:
%   - RMS uncertainty across all fits (sigma_a_total, sigma_b_total)
%   - Individual uncertainties (da_all, db_all) for each [Jcap1, Jcap2] pair
%   - Jcap_min, Jcap_max: min/max(Jcap1, Jcap2) for each fit
%   - ab_min, ab_max: min/max(a, b) for each fit
%
% NaN uncertainties are replaced with 1e-3 to keep all points for plots.
% -------------------------------------------------------------------------

% === NEW: pull uncertainties directly from the MAT table =============
T = load(matFilePath).fitTable;     % matFilePath is the first argument

da_all = mean([T.daPlus , T.daMinus],2,'omitnan');   % symmetric σ
db_all = mean([T.dbPlus , T.dbMinus],2,'omitnan');

a_all  = T.a;
b_all  = T.b;

Jcap_min = min([T.Jcap1 , T.Jcap2],[],2);
Jcap_max = max([T.Jcap1 , T.Jcap2],[],2);

ab_min = min([a_all , b_all],[],2);
ab_max = max([a_all , b_all],[],2);
% =====================================================================


    % --- Replace NaNs in uncertainty with fallback value ---
    fallback = 1e-3;
    nan_a = isnan(da_all);
    nan_b = isnan(db_all);
    if any(nan_a) || any(nan_b)
        fprintf('⚠ Replacing %d NaN da and %d NaN db with fallback = %.1e\n', ...
            sum(nan_a), sum(nan_b), fallback);
    end
    da_all(nan_a) = fallback;
    db_all(nan_b) = fallback;

    % No filtering – keep all entries (even those that had NaN before)
    
    % RMS propagated uncertainty
    sigma_a_total = sqrt(mean(da_all.^2));
    sigma_b_total = sqrt(mean(db_all.^2));

    % Compute min/max(a,b) for each pair
    ab_min = min([a_all, b_all], [], 2);
    ab_max = max([a_all, b_all], [], 2);

    % Log results
    fprintf('\nOverall RMS uncertainty:\n');
    fprintf('  σ_a = %.3e\n', sigma_a_total);
    fprintf('  σ_b = %.3e\n', sigma_b_total);

    % Create struct
    uncertStruct = struct();
    uncertStruct.sigma_a_total = sigma_a_total;
    uncertStruct.sigma_b_total = sigma_b_total;
    uncertStruct.da_all = da_all;
    uncertStruct.db_all = db_all;
    uncertStruct.a_all = a_all;
    uncertStruct.b_all = b_all;
    uncertStruct.Jcap_min = Jcap_min;
    uncertStruct.Jcap_max = Jcap_max;
    uncertStruct.ab_min = ab_min;
    uncertStruct.ab_max = ab_max;

    % Save to MAT file
    save(fullfile(output_path, 'overall_uncertainty_data.mat'), 'uncertStruct');
end

