close all
clear all 
clc

output_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon_cylinder\V2/my_post_proc_figures\Keff\no_coefficiants\';
mat_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon_cylinder\V2\my_post_proc_figures\Keff\no_coefficiants\postproc Keff fit results.mat';
postStructmat_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon_cylinder\V2\my_post_proc_figures\postProcessingStruct_Intrinsic.mat';

Jcap1 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
Jcap2 = [0.15, 0.135, 0.12, 0.105, 0.09, 0.076, 0.061, 0.046, 0.031, 0.016, 0.005, 0.001];
filter_AB = false;
min_PSI_points = 1;
R2_cutoff = 0.001;

run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB, R2_cutoff)
plot_fit_results_from_TXT(mat_path);

function run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB, R2_cutoff)
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
                postproc_fit_results(N1, N2, postStructmat_path, output_path, filter_AB, min_PSI_points, R2_cutoff);
            catch ME
                warning('Failed processing pair [%.4f, %.4f]: %s', N1, N2, ME.message);
            end
        end
    end
end



function postproc_fit_results(Jcap1, Jcap2, postStructmat_path, ...
                              output_path, filter_AB, ...
                              min_PSI_points, R2_cutoff)
% POSTPROC_FIT_RESULTS
% 1) Takes Jcap1, Jcap2, and a .mat struct file and extracts data directly
%    corresponding to the mixed-cap pair [Jcap1, Jcap2].
% 2) Retrieves [Jcap1, Jcap1] and [Jcap2, Jcap2] if available.
% 3) Performs fitting with model: Keff = A / (1 - phi)
% 4) Appends the result to postproc Keff fit results.mat (struct array)
%    and saves a .fig of the fit.

% -------------------------------------------------------------------------
% initialise
% -------------------------------------------------------------------------
if ~exist(output_path, 'dir');  mkdir(output_path);  end
S = load(postStructmat_path);
fn = fieldnames(S);
if numel(fn) ~= 1 || ~isstruct(S.(fn{1}))
    error('File %s must contain ONE struct variable.', postStructmat_path);
end
S = S.(fn{1});

tol = 1e-6;   aVal = Jcap1;   bVal = Jcap2;
if abs(aVal - bVal) < tol
    error('Jcap1 and Jcap2 must be different.');
end

% -------------------------------------------------------------------------
% gather AB / AA / BB data
% -------------------------------------------------------------------------
phiAB = []; psiAB = []; Z_AB = [];
dataAA = struct('phi',[],'psi',[],'Z',[]);
dataBB = dataAA;

for k = 1:numel(S)
    jcaps = unique(S(k).Jcap(:));   jcaps = sort(jcaps);
    switch numel(jcaps)
        case 2  % AB
            if all(ismembertol([aVal bVal],jcaps,tol))
                phiAB = [phiAB; S(k).phi(:)];
                psiAB = [psiAB; S(k).psi(:)];
                Z_AB  = [Z_AB ; S(k).keff(:)];
            end
        case 1  % AA or BB
            if abs(jcaps - aVal) < tol
                dataAA.phi = [dataAA.phi; S(k).phi(:)];
                dataAA.Z   = [dataAA.Z  ; S(k).keff(:)];
            elseif abs(jcaps - bVal) < tol
                dataBB.phi = [dataBB.phi; S(k).phi(:)];
                dataBB.Z   = [dataBB.Z  ; S(k).keff(:)];
            end
    end
end
if isempty(phiAB)
    error('Pair [%.4f, %.4f] not found in %s', aVal, bVal, postStructmat_path);
end

% raw copies (for plotting later)
phiAB_raw = phiAB;  psiAB_raw = psiAB;  Z_AB_raw = Z_AB;

% -------------------------------------------------------------------------
% optional AB filtering
% -------------------------------------------------------------------------
phiAB_filt = phiAB;  psiAB_filt = psiAB;  Z_AB_filt = Z_AB;
if filter_AB
    [uniquePsi, ~, idPsi] = uniquetol(psiAB, 1e-6);
    psiCounts  = accumarray(idPsi, 1);
    goodPsi    = psiCounts(idPsi) >= min_PSI_points;
    phiAB = phiAB(goodPsi);  psiAB = psiAB(goodPsi);  Z_AB = Z_AB(goodPsi);

    tolPsi = 1e-6;  minRSq = 0.9;  minPts = 10;
    [phiAB_filt, psiAB_filt, Z_AB_filt] = ...
        filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, ...
                                     tolPsi, minRSq, minPts);
end

% -------------------------------------------------------------------------
% build combined dataset (psi forced symmetric)
% -------------------------------------------------------------------------
if aVal < bVal
    dataAA.psi = ones(size(dataAA.phi));
    dataBB.psi = zeros(size(dataBB.phi));
else
    dataAA.psi = zeros(size(dataAA.phi));
    dataBB.psi = ones(size(dataBB.phi));
end
psiAll = [psiAB_filt; dataAA.psi; dataBB.psi];
phiAll = [phiAB_filt; dataAA.phi; dataBB.phi];
ZAll   = [Z_AB_filt ; dataAA.Z   ; dataBB.Z ];
if isempty(psiAll);  return;  end

% -------------------------------------------------------------------------
% fit
% -------------------------------------------------------------------------
[r2, A_fit, SE_A, b_fit, SE_b] = fitWithFixedKeffModel(psiAll, phiAll, ZAll);

% -------------------------------------------------------------------------
% write / append to .mat   (field-sync patch)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% simple file handling  – erase once, then just append
% -------------------------------------------------------------------------
persistent firstCall
if isempty(firstCall);  firstCall = true;  end

outMat = fullfile(output_path, 'postproc Keff fit results.mat');

if firstCall                     % happens only on the very first call
    if exist(outMat,'file');  delete(outMat);  end     % start fresh
    fitResults = struct([]);                        % empty struct array
    firstCall  = false;
else                            % subsequent calls in the same MATLAB run
    if exist(outMat,'file')
        tmp = load(outMat,'fitResults');
        fitResults = tmp.fitResults;
    else
        fitResults = struct([]);
    end
end


thisRes = struct( ...
    'pairAB'      , [aVal bVal], ...
    'pairAA'      , [aVal aVal], ...
    'pairBB'      , [bVal bVal], ...
    'A'           , A_fit, ...
    'SE_A'        , SE_A, ...
    'b'           , b_fit, ...
    'SE_b'        , SE_b, ...
    'R2'          , r2, ...
    'fitEquation' , {'Keff = 1/(1 - \phi)'}, ...
    'numRawAB'    , numel(Z_AB_raw), ...
    'numKeptAB'   , numel(Z_AB_filt), ...
    'timestamp'   , datetime('now') ...
);




% ---- make sure *all* structs have the same fields ----------------------
allFlds = union(fieldnames(fitResults), fieldnames(thisRes));
for f = 1:numel(allFlds)
    fld = allFlds{f};
    if ~isfield(fitResults, fld)
        % add missing field to every element, keep shape:
        [fitResults.(fld)] = deal([]);
    end
    if ~isfield(thisRes, fld)
        % add missing field to thisRes
        thisRes.(fld) = [];
    end
end

% ---- now safe to concatenate ------------------------------------------
fitResults = [fitResults; thisRes];      %#ok<AGROW>
save(outMat, 'fitResults');              % overwrite / create
fprintf('Appended result to %s  (total fits: %d)\n', outMat, numel(fitResults));

% -------------------------------------------------------------------------
% figure (unchanged, except legend uses cleaned labels)
% -------------------------------------------------------------------------
hFig = figure('Visible','on');  hold on;  grid on;  view(3);
scatter3(phiAB_raw, psiAB_raw, Z_AB_raw, 36,'d', ...
         'MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','none', ...
         'DisplayName',sprintf('[%.4f,%.4f] AB – raw',aVal,bVal));
scatter3(phiAB_filt,psiAB_filt,Z_AB_filt,40,'d', ...
         'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor','none', ...
         'DisplayName',sprintf('[%.4f,%.4f] AB – filtered',aVal,bVal));
scatter3(phiAll,psiAll,ZAll,45,'d','MarkerEdgeColor','k','MarkerFaceColor','k', ...
         'DisplayName','All inliers');
if ~isempty(dataAA.psi)
    scatter3(dataAA.phi,dataAA.psi,dataAA.Z,60,'^','MarkerEdgeColor','r', ...
             'MarkerFaceColor','r','DisplayName',sprintf('[%.4f,%.4f] AA',aVal,aVal));
end
if ~isempty(dataBB.psi)
    scatter3(dataBB.phi,dataBB.psi,dataBB.Z,60,'s','MarkerEdgeColor','b', ...
             'MarkerFaceColor','b','DisplayName',sprintf('[%.4f,%.4f] BB',bVal,bVal));
end
phi_lin = linspace(min(phiAll), max(phiAll), 50);
psi_lin = linspace(min(psiAll), max(psiAll), 50);
[phi_grid, psi_grid] = meshgrid(phi_lin, psi_lin);

Z_grid_fit = 1 ./ (1 - phi_grid);
Z_grid_fit(phi_grid >= 1) = NaN;

surf(phi_grid, psi_grid, Z_grid_fit, 'EdgeColor','none', ...
     'FaceAlpha',0.5,'FaceColor','interp', ...
     'DisplayName','K_{eff} theory = 1/(1-\phi)');



xlabel('\phi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
ylabel('\psi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
zlabel('K_{eff} [K_{B}T]',   'FontSize',40, 'Color','red', 'FontWeight','bold');
title(sprintf('K_{eff} fit [J_a,J_b]=[%.4f,%.4f]', aVal,bVal),   'FontSize',30, 'Color','red', 'FontWeight','bold');
subtitle('Keff - spherical model');
legend('Location','best');  hold off;
savefig(hFig, fullfile(output_path, ...
        sprintf('fit_phi_psi_%.4f_%.4f.fig',aVal,bVal)));
close(hFig);
end



function plot_fit_results_from_TXT(matFilePath)
% PLOT_FIT_RESULTS_FROM_TXT  – Plot R² and A vs mismatch Δ
% Works with 'postproc Keff fit results.mat' containing variable fitResults.
%
% Example:
%   plot_fit_results_from_TXT('D:\...\postproc Keff fit results.mat')

    % ---------------- load ------------------------------------------------
    data = load(matFilePath, 'fitResults');
    if ~isfield(data, 'fitResults') || isempty(data.fitResults)
        error('No ''fitResults'' variable found in %s', matFilePath);
    end
    R = data.fitResults;

    % ---------------- collect Δ, R², A -----------------------------------
    n = numel(R);
    deltaVals = zeros(1, n);
    r2Vals    = zeros(1, n);
    labels    = cell(1, n);


    for k = 1:n
        lowVal   = min(R(k).pairAA(1), R(k).pairBB(1));
        highVal  = max(R(k).pairAA(1), R(k).pairBB(1));
        deltaVals(k) = highVal - lowVal;
        r2Vals(k) = R(k).R2;
        labels{k} = sprintf('[%.4g, %.4g]', lowVal, highVal);
    end


    % ---------------- sort for tidy plots ---------------------------------
    [deltaVals, idx] = sort(deltaVals);
    r2Vals = r2Vals(idx);
    labels = labels(idx);

    [folderPath, ~, ~] = fileparts(matFilePath);

    % ---------------- R² vs Δ figure --------------------------------------
    hFig = figure('Visible', 'on');
    plot(deltaVals, r2Vals, 'cd-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Mismatch \Delta (High - Low)',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    ylabel('R^2',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    title('R^2 vs. Mismatch \Delta',   'FontSize',30, 'Color','red', 'FontWeight','bold');
    grid on;
    for i = 1:numel(deltaVals)
        text(deltaVals(i), r2Vals(i), labels{i}, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right', ...
            'FontSize', 8, 'Color', 'k');
    end
    savefig(hFig, fullfile(folderPath, 'R2Values_vs_Delta.fig'));
    close(hFig);

end



function [R_squared, A_fit, SE_A, b_fit, SE_b] = fitWithFixedKeffModel(psiAll, phiAll, ZAll)
    % Fixed model: Keff = 1/(1 - phi)
    % We only evaluate goodness-of-fit (R^2) to this surface.
    % For pipeline compatibility we return A=1, b=1, SEs=0.

    % Safety
    ZAll(ZAll <= 0)        = eps;
    phiAll(phiAll <= 0)    = eps;
    phiAll(phiAll >= 0.9999) = 0.9999;

    % Predicted Keff from the fixed model
    Zhat = 1 ./ (1 - phiAll(:));

    % R^2 in original Z space
    SS_res = sum((ZAll(:) - Zhat).^2);
    SS_tot = sum((ZAll(:) - mean(ZAll(:))).^2);
    R_squared = 1 - SS_res/SS_tot;

    % Fixed parameters for compatibility
    A_fit = 1;
    b_fit = 1;
    SE_A  = 0;
    SE_b  = 0;

    fprintf('Fixed model  Keff = 1/(1-phi):  R² = %.4f\n', R_squared);
end





function [phiAB_filt, psiAB_filt, Z_AB_filt] = ...
    filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, ...
                                  tolPsi, minRSquare, minPointsPerPsi)

    if nargin < 4, tolPsi = 1e-6; end
    if nargin < 5, minRSquare = 0.9; end
    if nargin < 6, minPointsPerPsi = 5; end

    residualThreshold = 0.99;
    maxAllowedOutlierFrac = 0.1;

    [uniquePsi, ~, idPsi] = uniquetol(psiAB, tolPsi);
    nBands = numel(uniquePsi);

    keepMask = false(size(psiAB));  % this will be returned directly

    for k = 1:nBands
        idx = (abs(psiAB - uniquePsi(k)) < tolPsi);
        if sum(idx) < minPointsPerPsi
            continue;
        end

        phi_k = phiAB(idx);
        Z_k   = Z_AB(idx);

        valid = isfinite(phi_k) & isfinite(Z_k);
        phi_k = phi_k(valid);
        Z_k   = Z_k(valid);

        if numel(Z_k) < minPointsPerPsi
            continue;
        end

        % Fit and residuals
        p = polyfit(phi_k, Z_k, 1);
        Z_fit = polyval(p, phi_k);
        residuals = abs(Z_k - Z_fit);

        outlierFraction = sum(residuals > residualThreshold) / numel(Z_k);

        SS_res = sum((Z_k - Z_fit).^2);
        SS_tot = sum((Z_k - mean(Z_k)).^2);
        R2 = 1 - SS_res / SS_tot;

        if outlierFraction <= maxAllowedOutlierFrac
            keepMask(idx) = true;  % mark those original indices
            status = 'KEPT';
        else
            status = 'DISCARDED';
        end

        fprintf('PSI = %.4f : R² = %.4f  |  Outlier fraction = %.2f → %s\n', ...
                uniquePsi(k), R2, outlierFraction, status);
    end

    % Apply final filtering
    phiAB_filt = phiAB(keepMask);
    psiAB_filt = psiAB(keepMask);
    Z_AB_filt  = Z_AB(keepMask);
    fprintf('→ Final kept points: %d out of %d\n', nnz(keepMask), numel(keepMask));
end

