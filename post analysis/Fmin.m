close all
clear all 
clc


output_path = 'C:\Users\User\OneDrive/my_output_directory/new_results_Fmin/cylinder/';
mat_path = 'C:\Users\User\OneDrive/my_output_directory/new_results_Fmin/cylinder\postproc Fmin PHI fit results.mat';
postStructmat_path = 'C:\Users\User\OneDrive/scripts\postProcessingStruct_Intrinsic_cylinder.mat';  % <-- Your custom folder path

Jcap1 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
Jcap2 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
min_PSI_points = 1;
filter_AB = false;
filter_AA_BB = false;
R2_cutoff_AA_BB = 0.2;
R2_cutoff = 0.7;


% run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB, filter_AA_BB, R2_cutoff, R2_cutoff_AA_BB);
plot_fit_results_from_MAT(mat_path);
% build_deltaJcap_scans(output_path)
% plot_deltaJcap_fit_results(output_path);


function run_for_all_Jcap_combinations(Jcap1, Jcap2, postStructmat_path, output_path, min_PSI_points, filter_AB, filter_AA_BB, R2_cutoff, R2_cutoff_AA_BB)
        if ~exist(output_path, 'dir')
            mkdir(output_path);
        end
% Ensure Jcap1 and Jcap2 are not empty
    if isempty(Jcap1) || isempty(Jcap2)
        error('Jcap1 and Jcap2 must not be empty');
    end

    processedPairs = containers.Map();  % keys like '0.001_0.150'

    % Loop through all combinations of Jcap1 and Jcap2
    for j = 1:length(Jcap1)
        for k = 1:length(Jcap2)
            N1 = Jcap1(j);
            N2 = Jcap2(k);

            min_cord = 25;
            large_Jcap = max(N1, N2);
            small_Jcap = min(N1, N2);
            delta_Jcap = large_Jcap-small_Jcap;

            if (delta_Jcap <= 0.01) && (large_Jcap > 0.04) && (small_Jcap > 0.04)
                min_cord = 10;
            end

            % Ensure they are not the same
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

            % Call your post-processing fitting function
            fprintf('Processing pair [%.4f, %.4f]\n', N1, N2);
            try
                postproc_fit_selected_pair_struct(postStructmat_path, output_path, N1, N2, min_PSI_points, filter_AB, filter_AA_BB, R2_cutoff, R2_cutoff_AA_BB);
            catch ME
                warning('Failed processing pair [%.4f, %.4f]: %s', N1, N2, ME.message);
            end
        end
    end
end


function postproc_fit_selected_pair_struct(structPath, output_path, Jcap1, Jcap2, min_PSI_points, filter_AB, filter_AA_BB, R2_cutoff, R2_cutoff_AA_BB)
% POSTPROC_FIT_SELECTED_PAIR_STRUCT
% -------------------------------------------------------------------------
% Build a 3-D plot (φ, ψ, Fmin) and perform a surface–fit for ONE specific
% mixed-cap pair [Jcap1 , Jcap2] **directly from the post-processing struct**.
%
% • AB  = mixed-cap   [Jcap1 , Jcap2]    (ψ varies – taken “as is”)
% • AA  = symmetric   [Jcap1 , Jcap1]    (ψ hard-coded)
% • BB  = symmetric   [Jcap2 , Jcap2]    (ψ hard-coded)
%
% The struct is assumed to contain (at minimum) the fields:
%     Jcap, phi, psi, Fmin, proteinNumber
% produced by the normal post-processing pipeline.
%
% Results are appended to   “…/postproc fit results.txt”
% and the figure is saved next to the *.mat file.
% -------------------------------------------------------------------------

tol   = 1e-6;
aVal  = Jcap1;
bVal  = Jcap2;
if abs(aVal-bVal)<tol
    error('Jcap1 and Jcap2 must be different.');
end

% --- Orientation tag and styling: AB vs BA -------------------------------
% Convention:
%   AB orientation : Jcap1 < Jcap2
%   BA orientation : Jcap1 > Jcap2
if aVal < bVal
    orientStr     = 'AB';
    colorAB_raw   = [0.0 0.0 0.0];      % black-ish for AB raw points
    markerAB_raw  = 'd';                % diamond markers for AB
else
    orientStr     = 'BA';
    colorAB_raw   = [0.0 0.45 0.74];    % bluish for BA raw points
    markerAB_raw  = 'o';                % circles for BA
end


%% ------------------------------------------------------------------------
%% 1) -------- load the struct --------------------------------------------
%% ------------------------------------------------------------------------
S = load(structPath);
fn = fieldnames(S);
if numel(fn)~=1 || ~isstruct(S.(fn{1}))
    error('File %s must contain ONE struct variable.',structPath);
end
S = S.(fn{1});        % the 1×N struct array

%% ------------------------------------------------------------------------
%% 2) -------- split into AB / AA / BB  (AB+BA mixed, orientation kept) ---
%% ------------------------------------------------------------------------
phiAB = [];  psiAB = [];  Z_AB  = [];  protAB = [];  jinAB = [];  kAB = [];
orientAB = [];   % 1 = AB orientation [aVal,bVal], 2 = BA orientation [bVal,aVal]

dataAA = struct('phi',[],'psi',[],'Z',[],'prot',[],'cord',[],'jin',[],'k',[]);
dataBB = struct('phi',[],'psi',[],'Z',[],'prot',[],'cord',[],'jin',[],'k',[]);

for k = 1:numel(S)
    % Jcap is stored as a 1x2 vector, order carries orientation
    pair = S(k).Jcap(:).';      % force row [J1 J2]
    if numel(pair) ~= 2
        warning('Entry %d has unexpected Jcap size, skipping.',k);
        continue;
    end
    j1 = pair(1);
    j2 = pair(2);

    % ----- mixed pair: either [aVal,bVal] or [bVal,aVal] -----------------
    sortedPair = sort(pair);
    if numel(unique(sortedPair))==2 && ...
       abs(sortedPair(1)-aVal)<tol && abs(sortedPair(2)-bVal)<tol

        nPts = numel(S(k).phi(:));

        phiAB  = [phiAB ; S(k).phi(:)];
        psiAB  = [psiAB ; S(k).psi(:)];
        Z_AB   = [Z_AB  ; S(k).Fmin(:)];
        protAB = [protAB; S(k).proteinNumber(:)];
        jinAB  = [jinAB ; S(k).Jintrin(:)];
        kAB    = [kAB   ; S(k).keff(:)];

        % orientation flag for each point
        if abs(j1-aVal)<tol && abs(j2-bVal)<tol
            thisOrient = 1;     % AB
        elseif abs(j1-bVal)<tol && abs(j2-aVal)<tol
            thisOrient = 2;     % BA
        else
            thisOrient = 0;     % should not happen, but be safe
        end
        orientAB = [orientAB; repmat(thisOrient, nPts, 1)];
        continue
    end

    % ----- AA: [aVal,aVal] ----------------------------------------------
    if abs(j1-aVal)<tol && abs(j2-aVal)<tol
        dataAA.phi  = [dataAA.phi ; S(k).phi(:)];
        dataAA.Z    = [dataAA.Z   ; S(k).Fmin(:)];
        dataAA.prot = [dataAA.prot; S(k).proteinNumber(:)];
        dataAA.jin  = [dataAA.jin ; S(k).Jintrin(:)];
        dataAA.k    = [dataAA.k   ; S(k).keff(:)];
        continue
    end

    % ----- BB: [bVal,bVal] ----------------------------------------------
    if abs(j1-bVal)<tol && abs(j2-bVal)<tol
        dataBB.phi  = [dataBB.phi ; S(k).phi(:)];
        dataBB.Z    = [dataBB.Z   ; S(k).Fmin(:)];
        dataBB.prot = [dataBB.prot; S(k).proteinNumber(:)];
        dataBB.jin  = [dataBB.jin ; S(k).Jintrin(:)];
        dataBB.k    = [dataBB.k   ; S(k).keff(:)];
        continue
    end
end

%% ------------------ AB filtering

% save the raw AB data
phiAB_raw  = phiAB;
psiAB_raw  = psiAB;
Z_AB_raw   = Z_AB;
protAB_raw = protAB;
jinAB_raw  = jinAB;     % NEW
kAB_raw    = kAB;       % NEW
orientAB_raw  = orientAB;   % new

% save the raw AA data
phiAA_raw  = dataAA.phi;
Z_AA_raw   = dataAA.Z;
protAA_raw = dataAA.prot;
jinAA_raw  = dataAA.jin;   % NEW
kAA_raw    = dataAA.k;     % NEW

% save the raw BB data
phiBB_raw  = dataBB.phi;
Z_BB_raw   = dataBB.Z;
protBB_raw = dataBB.prot;
jinBB_raw  = dataBB.jin;   % NEW
kBB_raw    = dataBB.k;     % NEW

if (filter_AB)

    % Filter AB
    phiAB = phiAB(maskAB);
    psiAB = psiAB(maskAB);
    Z_AB  = Z_AB(maskAB);
    protAB = protAB(maskAB);

    %% ----------- Filter AB data by ψ band population -------------------------
    tolPsi       = 1e-6;   % tolerance to group ψ values
    minPointsPerPsi = min_PSI_points;   % keep ψ bands with at least 20 points
    
    [uniquePsi,~,idPsi] = uniquetol(psiAB, tolPsi);
    psiCounts = accumarray(idPsi, 1);
    
    % Logical mask: keep only points where the ψ band has enough members
    validBands = psiCounts(idPsi) >= minPointsPerPsi;
    
    % Apply mask to AB
    phiAB  = phiAB(validBands);
    psiAB  = psiAB(validBands);
    Z_AB   = Z_AB(validBands);
    protAB = protAB(validBands);
    orientAB = orientAB(maskAB);      % NEW
    midJcapAB = midJcapAB(validBands);
    cordAB = cordAB(validBands);
    jinAB  = jinAB(validBands);
    kAB    = kAB(validBands);
    
    [phiAB, psiAB, Z_AB, protAB] = filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, protAB, 1e-6, R2_cutoff, 5);
    if isempty(phiAB)
        error('AB data is empty after filtering (min_cord, min_PSI_points, R² check). Check your thresholds or input data.');
    end
    kept = [phiAB, psiAB, Z_AB, protAB];
    raw  = [phiAB_raw, psiAB_raw, Z_AB_raw, protAB_raw];
    [~,loc] = ismember(kept, raw, 'rows');
    jinAB = jinAB_raw(loc);
    kAB   = kAB_raw(loc);

end


if (filter_AA_BB)
    %% ------------------------------------------------------------------
    % Filter AA
    maskAA = dataAA.cord >= min_cord;
    dataAA.phi = dataAA.phi(maskAA);
    dataAA.Z   = dataAA.Z(maskAA);
    dataAA.prot = dataAA.prot(maskAA);
    dataAA.jin = dataAA.jin(maskAA);
    dataAA.k   = dataAA.k(maskAA);
    % Rebuild psi as done before
    dataAA.psi = zeros(size(dataAA.phi));  

    [dataAA.phi, dataAA.Z, dataAA.prot] = ...
        filter_by_protein_fit_quality(dataAA.phi, dataAA.Z, dataAA.prot, R2_cutoff_AA_BB, 1);
    
    % Filter BB
    maskBB = dataBB.cord >= min_cord;
    dataBB.phi = dataBB.phi(maskBB);
    dataBB.Z   = dataBB.Z(maskBB);
    dataBB.prot = dataBB.prot(maskBB);
    dataBB.psi = ones(size(dataBB.phi));
    dataBB.jin = dataBB.jin(maskBB);
    dataBB.k   = dataBB.k(maskBB);

    [dataBB.phi, dataBB.Z, dataBB.prot] = ...
        filter_by_protein_fit_quality(dataBB.phi, dataBB.Z, dataBB.prot, R2_cutoff_AA_BB, 1);
end

if isempty(phiAB)
    error('Pair [%.4f, %.4f] not found in %s',aVal,bVal,structPath);
end

% hard-code ψ for symmetric sets
if aVal < bVal
    dataAA.psi = ones(size(dataAA.phi));   % low Jcap → ψ = 0
    dataBB.psi = zeros (size(dataBB.phi));   % high     → ψ = 1
else
    dataAA.psi = zeros (size(dataAA.phi));
    dataBB.psi = ones (size(dataBB.phi));
end


%% ------------------------------------------------------------------------
%% 3) -------- fit on AB inliers ------------------------------------------
%% ------------------------------------------------------------------------

[fitRes,gof,psiC,phiC,ZC,idxIn,useLog,Zmax] = ...
    fitWithMultipleStartPoints( ...
        aVal,bVal, psiAB,phiAB,Z_AB,   ...   % AB (first!)
        dataAA.psi,dataAA.phi,dataAA.Z,              ...   % AA  (optional)
        dataBB.psi,dataBB.phi,dataBB.Z);                  % BB  (optional)

phiIn = phiAB(idxIn);
psiIn = psiAB(idxIn);
ZIn   = Z_AB(idxIn);
orientIn = orientAB(idxIn);           % orientation for inliers

% --- cache AB inliers per pair for cross-ΔJcap plots (self-contained) ---
try
    % order-preserving tag (use sort([aVal bVal]) if you want canonical tag)
    pairTag    = sprintf('%.4f_%.4f', aVal, bVal);

    % ensure per-pair folder exists *now*, independent of later code
    if ~exist(output_path,'dir'), mkdir(output_path); end
    pairFolder = fullfile(output_path, ['pair_' pairTag]);
    if ~exist(pairFolder,'dir'), mkdir(pairFolder); end

    % package compact data
    deltaJcap = abs(aVal - bVal);
    pairData = struct( ...
        'Jcap1', aVal, ...
        'Jcap2', bVal, ...
        'deltaJcap', deltaJcap, ...
        'phi',   phiIn(:), ...
        'psi',   psiIn(:), ...
        'Fmin',  ZIn(:) );

    % save cache
    pairDataFile = fullfile(pairFolder, sprintf('pair_AB_inliers_%s.mat', pairTag));
    save(pairDataFile, 'pairData', '-v7');

catch ME
    warning('AB-inliers cache skipped for [%g,%g]: %s', aVal, bVal, ME.message);
end



% Axis-limits helpers
psiAll = [psiAB ; dataAA.psi ; dataBB.psi];
phiAll = [phiAB ; dataAA.phi ; dataBB.phi];

%% 4) -------- append results  (same text, but into a MAT file) ------------

[outFolder,~,~] = fileparts(output_path);
if ~exist(outFolder,'dir'), mkdir(outFolder); end

matFile = fullfile(outFolder,'postproc Fmin PHI fit results.mat');

% ---- per-pair subfolder (order-preserving; change to sort([aVal bVal]) if you prefer) ----
pairTag    = sprintf('%.4f_%.4f', aVal, bVal);
pairFolder = fullfile(outFolder, ['pair_' pairTag]);
if ~exist(pairFolder,'dir'), mkdir(pairFolder); end


% ---------- load or initialise ------------------------------------------
equationStr = "Fit Equation: Fmin = (Kmem/2)*((psi*(1-psi)))*(phi^a)*b";

if isfile(matFile)
    S = load(matFile);                % may contain fitTable & equationStr
    if isfield(S,'fitTable')
        fitTable = S.fitTable;
    else
fitTable = table([],[],[],[],[],[], ...
    'VariableNames',{'Jcap1','Jcap2','a','b','Rsq','a_err'});
    end
else
fitTable = table([],[],[],[],[],[], ...
    'VariableNames',{'Jcap1','Jcap2','a','b','Rsq','a_err'});
end

% ---------- append numeric row ------------------------------------------
% ---- alpha error from 95% CI ----
try
    CI = confint(fitRes);        % [a_low b_low ; a_high b_high]
    a_err = abs(CI(2,1) - CI(1,1))/2;
catch
    a_err = NaN;
end

newRow = {aVal , bVal , fitRes.a , fitRes.b , gof.rsquare , a_err};
fitTable = [fitTable ; newRow];


% ---------- save (add equationStr once) ---------------------------------
save(matFile,'fitTable','equationStr','-v7');        % use -v7.3 if >2 GB

fprintf('Fit results appended to MAT-file:\n   %s\n', matFile);


%% ------------------------------------------------------------------------
%% 5) -------- colours for proteins ---------------------------------------
%% ------------------------------------------------------------------------
protSet = unique([dataAA.prot ; dataBB.prot]);
nProt = max(protSet);               % how many proteins?
hues = linspace(0,1,nProt+1).';     % drop the last to avoid red-wrap
hues(end) = [];
cMap = hsv2rgb([hues , 0.75*ones(nProt,1) , 0.9*ones(nProt,1)]);

%% ------------------------------------------------------------------------
%% 6) -------- plotting ----------------------------------------------------
%% ------------------------------------------------------------------------
hFig = figure('Visible','on');
set(hFig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on;  grid on;  view(3);
set(findall(hFig,'Type','axes'),'FontSize',28,'LineWidth',1.2);

% AB and BA raw, marked separately
maskAB_orient = (orientAB_raw == 1);
maskBA_orient = (orientAB_raw == 2);

% if any(maskAB_orient)
%     scatter3(phiAB_raw(maskAB_orient), ...
%              psiAB_raw(maskAB_orient), ...
%              Z_AB_raw(maskAB_orient), ...
%              40, [0 0 0], 'd', ...                     % black diamonds
%              'MarkerFaceColor','none','MarkerEdgeAlpha',0.3, ...
%              'DisplayName', sprintf('[%.4f, %.4f] (AB – all raw)', aVal, bVal));
% end
% 
% if any(maskBA_orient)
%     scatter3(phiAB_raw(maskBA_orient), ...
%              psiAB_raw(maskBA_orient), ...
%              Z_AB_raw(maskBA_orient), ...
%              40, [0 0.45 0.74], 'o', ...               % blue circles
%              'MarkerFaceColor','none','MarkerEdgeAlpha',0.3, ...
%              'DisplayName', sprintf('[%.4f, %.4f] (BA – all raw)', aVal, bVal));
% end


% --- OUTLIERS in AA ---
[~, idxKeepAA] = ismember([dataAA.phi, dataAA.Z, dataAA.prot], ...
                          [phiAA_raw, Z_AA_raw, protAA_raw], 'rows');
idxOutAA = setdiff(1:numel(phiAA_raw), idxKeepAA);

% Plot AA outliers as empty red diamonds
if ~isempty(idxOutAA)
    scatter3(phiAA_raw(idxOutAA), ...
             dataAA.psi(1)*ones(size(idxOutAA)), ...
             Z_AA_raw(idxOutAA), ...
             45, 'r', 'd', ...
             'MarkerFaceColor','none', 'MarkerEdgeAlpha',0.3, ...
             'DisplayName', sprintf('[%.4f, %.4f] (AA – outliers)', aVal, aVal));
end


% --- OUTLIERS in BB ---
[~, idxKeepBB] = ismember([dataBB.phi, dataBB.Z, dataBB.prot], ...
                          [phiBB_raw, Z_BB_raw, protBB_raw], 'rows');
idxOutBB = setdiff(1:numel(phiBB_raw), idxKeepBB);

% Plot BB outliers as empty blue diamonds
if ~isempty(idxOutBB)
    scatter3(phiBB_raw(idxOutBB), ...
             dataBB.psi(1)*ones(size(idxOutBB)), ...
             Z_BB_raw(idxOutBB), ...
             45, 'b', 'd', ...
             'MarkerFaceColor','none', 'MarkerEdgeAlpha',0.3, ...
             'DisplayName', sprintf('[%.4f, %.4f] (BB – outliers)', bVal, bVal));
end

% AB + BA – inliers, colour by protein, marker by orientation
protIn  = protAB(idxIn);
protSet = unique([dataAA.prot ; dataBB.prot ; protIn]);
cMap    = hsv(max(protSet));

assert( ...
    all( abs(psiIn - psiAB(idxIn))./max(psiAB) < 1e-12 ), ...
    'idxIn is still mis-aligned!' );

maskIn_AB = (orientIn == 1);
maskIn_BA = (orientIn == 2);

% AB inliers → same shape as AA ('^')
if any(maskIn_AB)
    addABInlierScatter( ...
        phiIn(maskIn_AB), psiIn(maskIn_AB), ZIn(maskIn_AB), ...
        protIn(maskIn_AB), ...
        '[%.4f, %.4f] (AB – inliers)', ...
        aVal, bVal, protSet, cMap, 1);   % 1 = AB → AA shape
end

% BA inliers → same shape as BB ('s')
if any(maskIn_BA)
    addABInlierScatter( ...
        phiIn(maskIn_BA), psiIn(maskIn_BA), ZIn(maskIn_BA), ...
        protIn(maskIn_BA), ...
        '[%.4f, %.4f] (BA – inliers)', ...
        aVal, bVal, protSet, cMap, 2);   % 2 = BA → BB shape
end



% AA & BB – colour-coded by protein
addSymmetricScatter(dataAA, '[%.4f, %.4f] (AA)', aVal, protSet, cMap, '^');
addSymmetricScatter(dataBB, '[%.4f, %.4f] (BB)', bVal, protSet, cMap, 's');

% fitted surface
psiGrid = linspace(min(psiAll), max(psiAll), 50);
phiGrid = linspace(min(phiAll), max(phiAll), 50);
[PsiMesh,PhiMesh] = meshgrid(psiGrid,phiGrid);
ZMesh = fitRes(PsiMesh,PhiMesh);

Zmax_surface = max(ZMesh(:));


surf(PhiMesh,PsiMesh,ZMesh, ...
     'EdgeColor','none','FaceAlpha',0.35, ...
     'DisplayName', sprintf('Fitted surface (%s)', orientStr));

% zlim([min(ZMesh(:)), Zmax_surface * 1.1]);
eqStr  = 'Fit Equation: F_{min} = (K_{mem}/2) (\psi(1-\psi)) \cdot (\phi^{\alpha}) \cdot \beta';
parStr = ['\alpha = ' num2str(fitRes.a,'%.3g') ', \beta = ' num2str(fitRes.b,'%.3g')];

xlabel('\phi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
ylabel('\psi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
zlabel('$F_{\min}\;\left[\frac{k_{\mathrm{B}}T}{\mathrm{nm}^{2}}\right]$', ...
           'Interpreter','latex', ...
           'FontSize',40, ...
           'Color','red', ...
           'FontWeight','bold');
title(sprintf('Fit for Jcap pair [%.4f, %.4f]', aVal, bVal), ...
      'FontSize', 20, 'Color', 'red', 'FontWeight', 'bold');
subtitle({eqStr, parStr}, ...
    'Interpreter','tex', 'FontSize',20, 'Color','black', 'FontWeight','bold');
legend('Location','best');
set(legend(),'FontSize',20);

savefig(hFig, fullfile(pairFolder, ...
        sprintf('Fmin_fit_%.4f_%.4f.fig',aVal,bVal)));
close(hFig);

% =====================  EXTRA FIGURE: Normalized Fmin  ====================
% Z_norm = Fmin / (k_eff * Jintrin^2)
safeDiv = @(Z,K,J) Z ./ (K .* (J.^2));

% Build normalized arrays with guards against bad values
Z_AB_norm = safeDiv(Z_AB, kAB, jinAB);
goodAB    = isfinite(Z_AB_norm) & (kAB>0) & (jinAB~=0);

% keep only the “good” AB points for plotting + sensitivity
phiAB_norm   = phiAB(goodAB);
psiAB_norm   = psiAB(goodAB);
ZAB_norm_good = Z_AB_norm(goodAB);
protAB_norm  = protAB(goodAB);
orientAB_norm = orientAB(goodAB);   % 1 = AB , 2 = BA

Z_AA_norm = safeDiv(dataAA.Z, dataAA.k, dataAA.jin);
goodAA    = isfinite(Z_AA_norm) & (dataAA.k>0) & (dataAA.jin~=0);

Z_BB_norm = safeDiv(dataBB.Z, dataBB.k, dataBB.jin);
goodBB    = isfinite(Z_BB_norm) & (dataBB.k>0) & (dataBB.jin~=0);


hFigN = figure('Visible','on');
set(hFigN, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
hold on; grid on; view(3);
set(findall(hFigN,'Type','axes'),'FontSize',28,'LineWidth',1.2);

% AB + BA normalized, shaped like AA/BB and coloured by protein
protSetNorm = unique(protAB_norm);
cMapNorm    = hsv(max(protSetNorm));   % or reuse cMap if you prefer

for p = protSetNorm(:).'
    % AB orientation (like AA): marker = '^'
    idxAB = (protAB_norm == p) & (orientAB_norm == 1);
    if any(idxAB)
        scatter3(phiAB_norm(idxAB), ...
                 psiAB_norm(idxAB), ...
                 ZAB_norm_good(idxAB), ...
                 40, cMapNorm(p,:), '^', 'filled', ...
                 'DisplayName', sprintf('[%.4f, %.4f] (AB – norm) proteins = %d', ...
                                         aVal, bVal, p));
    end

    % BA orientation (like BB): marker = 's'
    idxBA = (protAB_norm == p) & (orientAB_norm == 2);
    if any(idxBA)
        scatter3(phiAB_norm(idxBA), ...
                 psiAB_norm(idxBA), ...
                 ZAB_norm_good(idxBA), ...
                 40, cMapNorm(p,:), 's', 'filled', ...
                 'DisplayName', sprintf('[%.4f, %.4f] (BA – norm) proteins = %d', ...
                                         aVal, bVal, p));
    end
end


% AA – normalized (colour by protein)
if any(goodAA)
    for p = unique(dataAA.prot(goodAA)).'
        idx = goodAA & (dataAA.prot==p);
        scatter3(dataAA.phi(idx), dataAA.psi(idx), Z_AA_norm(idx), ...
                 40, [0.85 0.2 0.2], '^', 'filled', ...
                 'DisplayName', sprintf('[%.4f, %.4f] (AA – norm) proteins=%d', aVal, aVal, p));
    end
end

% BB – normalized (colour by protein)
if any(goodBB)
    for p = unique(dataBB.prot(goodBB)).'
        idx = goodBB & (dataBB.prot==p);
        scatter3(dataBB.phi(idx), dataBB.psi(idx), Z_BB_norm(idx), ...
                 40, [0.2 0.4 0.9], 's', 'filled', ...
                 'DisplayName', sprintf('[%.4f, %.4f] (BB – norm) proteins=%d', bVal, bVal, p));
    end
end

xlabel('\phi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
ylabel('\psi',   'FontSize',40, 'Color','red', 'FontWeight','bold');
zlabel('$\displaystyle \frac{F_{\min}}{k_{\mathrm{eff}}\,J_{\mathrm{int}}^{2}}$', ...
       'Interpreter','latex', ...
       'FontSize',30, 'Color','red', 'FontWeight','bold');
title(sprintf('Ratio of F_{min} (C) and theory constant k_{eff} J_{intrin}^2 for pair [%.4f, %.4f]', ...
      aVal, bVal), ...
      'FontSize',30, 'Color','red', 'FontWeight','bold');
legend('Location','best');
set(legend(),'FontSize',20);

savefig(hFigN, fullfile(pairFolder, ...
        sprintf('Fmin_NORMALIZED_fit_%.4f_%.4f.fig',aVal,bVal)));
close(hFigN);
% ==========================================================================

% =====================  QUANTIFY φ vs ψ sensitivity  ======================
% We’ll use AB (mixed) points only; AA/BB sit at ψ=1/0 and are edge cases.
phiABn = phiAB(goodAB);
psiABn = psiAB(goodAB);
ZABn   = Z_AB_norm(goodAB);

if ~isempty(phiABn)
    % --- group by ψ (within tolerance) ---
    tolPsi = 1e-6;
    [psiU,~,idPsi] = uniquetol(psiABn, tolPsi);
    
    % containers for per-ψ stats
    Srows = [];
    for kpsi = 1:numel(psiU)
        m = (idPsi == kpsi);
        phi_k = phiABn(m);
        Z_k   = ZABn(m);

        % sort by φ for stable plots
        [phi_k, ord] = sort(phi_k);
        Z_k = Z_k(ord);

        % linear fit Z ~ c0 + c1*φ
        if numel(phi_k)>=2 && std(phi_k)>0
            p1 = polyfit(phi_k, Z_k, 1);
            Zhat1 = polyval(p1, phi_k);
            R2_1 = 1 - sum((Z_k - Zhat1).^2) / max(eps, sum((Z_k - mean(Z_k)).^2));
            slope1 = p1(1);
        else
            R2_1 = NaN; slope1 = NaN;
        end

        % quadratic fit Z ~ q0 + q1*φ + q2*φ^2 (captures mild curvature in φ)
        if numel(phi_k)>=3 && std(phi_k)>0
            p2 = polyfit(phi_k, Z_k, 2);
            Zhat2 = polyval(p2, phi_k);
            R2_2 = 1 - sum((Z_k - Zhat2).^2) / max(eps, sum((Z_k - mean(Z_k)).^2));
        else
            R2_2 = NaN; p2 = [NaN NaN NaN];
        end

        % rank correlation (monotonicity) between Z and φ
        if numel(phi_k)>=3
            % robust Spearman (skip tiny/constant bands; fallback if toolbox missing)
            if numel(phi_k) >= 3 && std(phi_k) > 0 && std(Z_k) > 0
                try
                    rho = corr(phi_k, Z_k, 'type','Spearman','rows','complete');
                catch
                    % Fallback: Pearson on ranks (no Statistics Toolbox required)
                    rk1 = tiedrank_fallback(phi_k);
                    rk2 = tiedrank_fallback(Z_k);
                    C   = corrcoef(rk1, rk2);
                    rho = C(1,2);
                end
            else
                rho = NaN;
            end

        else
            rho = NaN;
        end

        % simple range-based sensitivity in φ within this ψ band
        dZ_dPhi_range = (max(Z_k)-min(Z_k)) / max(eps, (max(phi_k)-min(phi_k)));

        Srows = [Srows; struct( ...
            'psi', psiU(kpsi), ...
            'N', numel(phi_k), ...
            'slope_lin', slope1, ...
            'R2_lin', R2_1, ...
            'R2_quad', R2_2, ...
            'rho_spearman', rho, ...
            'dZ_dPhi_range', dZ_dPhi_range, ...
            'Z_mean', mean(Z_k), ...
            'Z_std',  std(Z_k) )];
    end

    % --- global variance decomposition: how much is φ vs ψ? ---------------
    % Treat ψ (band id) as categorical; φ as numeric.
    % Compute η²_ψ = between-ψ variance / total variance (ANOVA-style).
    Ztot = ZABn;
    mu   = mean(Ztot);
    SST  = sum((Ztot - mu).^2);

    % between-ψ
    SSpsi = 0;
    for kpsi = 1:numel(psiU)
        m = (idPsi == kpsi);
        nk = sum(m);
        if nk>0
            mk = mean(ZABn(m));
            SSpsi = SSpsi + nk * (mk - mu)^2;
        end
    end
    eta2_psi = SSpsi / max(eps, SST);

    % rough φ contribution after removing ψ means (within-band)
    % per-ψ means (length = numel(psiU)), then expand back to pointwise length
    idPsi  = double(idPsi(:));                                   % ensure column, double
    bandMean = accumarray(idPsi, ZABn, [numel(psiU) 1], @mean, NaN);
    Zcenter  = ZABn - bandMean(idPsi);

    goodC   = isfinite(Zcenter);
    if any(goodC)
        x = phiABn(goodC);  y = Zcenter(goodC);
        x = x(:); y = y(:);
        X = [ones(size(x)) x];          % intercept + slope
        beta = X \ y;                    % OLS
        yhat = X*beta;
        SSE = sum((y - yhat).^2);
        SST = sum((y - mean(y)).^2);
        R2_phi_given_psi = 1 - SSE/max(eps,SST);
    else
        R2_phi_given_psi = NaN;
    end


    % print a compact report
    fprintf('\n[φ vs ψ sensitivity | pair %.4f–%.4f]\n', aVal, bVal);
    fprintf('  η^2_ψ (variance explained by ψ bands)     : %.3f\n', eta2_psi);
    fprintf('  R^2(φ | ψ) (φ explains residual within ψ) : %.3f\n', R2_phi_given_psi);

    % table output to console (sorted by ψ)
    Tpsi = struct2table(Srows);
    Tpsi = sortrows(Tpsi, 'psi');
    disp(Tpsi);

%% ------------------------------------------------------------------------
%% EXPORTS: per-pair CSV and MAT + append to master CSV
%% ------------------------------------------------------------------------
pairTag = sprintf('%0.4f_%0.4f', aVal, bVal);

% -- make sure the per-pair folder exists
if ~exist(pairFolder,'dir'), mkdir(pairFolder); end

% --- Per-ψ table with pair identifiers
Tpsi_export = Tpsi;
Tpsi_export.pair_Jcap1 = repmat(aVal, height(Tpsi_export), 1);
Tpsi_export.pair_Jcap2 = repmat(bVal, height(Tpsi_export), 1);
Tpsi_export = movevars(Tpsi_export, {'pair_Jcap1','pair_Jcap2'}, 'Before', 1);

% --- Per-pair CSV (overwrites for this pair)
csvPair = fullfile(pairFolder, sprintf('Znorm_sensitivity_%s.csv', pairTag));
try
    writetable(Tpsi_export, csvPair);
catch
    fid = fopen(csvPair,'w');
    fprintf(fid, 'pair_Jcap1,pair_Jcap2,psi,N,slope_lin,R2_lin,R2_quad,rho_spearman,dZ_dPhi_range,Z_mean,Z_std\n');
    for r = 1:height(Tpsi_export)
        fprintf(fid, '%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
            Tpsi_export.pair_Jcap1(r), Tpsi_export.pair_Jcap2(r), Tpsi_export.psi(r), ...
            Tpsi_export.N(r), Tpsi_export.slope_lin(r), Tpsi_export.R2_lin(r), ...
            Tpsi_export.R2_quad(r), Tpsi_export.rho_spearman(r), ...
            Tpsi_export.dZ_dPhi_range(r), Tpsi_export.Z_mean(r), Tpsi_export.Z_std(r));
    end
    fclose(fid);
end

% --- Build a human-readable summary and store it inside the MAT (no .txt)
summaryLines = {};
summaryLines{end+1} = sprintf('[phi vs psi sensitivity | pair %0.4f–%0.4f]', aVal, bVal);
summaryLines{end+1} = sprintf('eta2_psi (variance explained by psi)      : %.6f', eta2_psi);
summaryLines{end+1} = sprintf('R2(phi | psi) (phi explains residual)     : %.6f', R2_phi_given_psi);
summaryLines{end+1} = '';
summaryLines{end+1} = 'Per-psi summary:';
summaryLines{end+1} = 'psi\tN\tslope_lin\tR2_lin\tR2_quad\trho_spear\tdZ/dphi_range\tZ_mean\tZ_std';
for r = 1:height(Tpsi)
    summaryLines{end+1} = sprintf('%.6g\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g', ...
        Tpsi.psi(r), Tpsi.N(r), Tpsi.slope_lin(r), Tpsi.R2_lin(r), ...
        Tpsi.R2_quad(r), Tpsi.rho_spearman(r), Tpsi.dZ_dPhi_range(r), ...
        Tpsi.Z_mean(r), Tpsi.Z_std(r));
end
summary_text = strjoin(summaryLines, newline);

% --- Compute collapse vector NOW (before saving & plotting)
Zcollapse = ZABn ./ bandMean(idPsi);
Zcollapse(~isfinite(Zcollapse)) = NaN;   % guard against 0/NaN means

% --- MAT dump (per pair; contains all intermediate arrays too)
matPair = fullfile(pairFolder, sprintf('Znorm_sensitivity_%s.mat', pairTag));
SensitivityReport = struct();
SensitivityReport.Jcap1 = aVal;
SensitivityReport.Jcap2 = bVal;
SensitivityReport.eta2_psi = eta2_psi;
SensitivityReport.R2_phi_given_psi = R2_phi_given_psi;
SensitivityReport.Tpsi = Tpsi;
SensitivityReport.psiU = psiU;
SensitivityReport.idPsi = idPsi;
SensitivityReport.phiAB = phiABn;
SensitivityReport.psiAB = psiABn;
SensitivityReport.ZnormAB = ZABn;
SensitivityReport.bandMean = bandMean;
SensitivityReport.Zcollapse = Zcollapse(:);
SensitivityReport.Zcenter   = Zcenter(:);
SensitivityReport.summary_text = summary_text;

try
    save(matPair, 'SensitivityReport', '-v7');
catch ME
    warning('Failed to save MAT file %s: %s', matPair, ME.message);
end

% --- Append to a master CSV (one row per ψ band per pair)
csvMaster = fullfile(outFolder, 'Znorm_sensitivity_MASTER.csv');
if ~isfile(csvMaster)
    writetable(Tpsi_export, csvMaster);   % create
else
    fid = fopen(csvMaster,'a');
    for r = 1:height(Tpsi_export)
        fprintf(fid, '%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n', ...
            Tpsi_export.pair_Jcap1(r), Tpsi_export.pair_Jcap2(r), Tpsi_export.psi(r), ...
            Tpsi_export.N(r), Tpsi_export.slope_lin(r), Tpsi_export.R2_lin(r), ...
            Tpsi_export.R2_quad(r), Tpsi_export.rho_spearman(r), ...
            Tpsi_export.dZ_dPhi_range(r), Tpsi_export.Z_mean(r), Tpsi_export.Z_std(r));
    end
    fclose(fid);
end

% ---- MASTER .MAT aggregator (append table rows across pairs)
matMaster = fullfile(outFolder, 'Znorm_sensitivity_MASTER.mat');
if isfile(matMaster)
    S = load(matMaster);            % expects S.masterTable if present
    if isfield(S,'masterTable')
        masterTable = S.masterTable;
    else
        masterTable = Tpsi_export([],:);  % empty same-schema table
    end
else
    masterTable = Tpsi_export([],:);      % initialize empty table
end

% Append current pair’s rows
masterTable = [masterTable; Tpsi_export];   %#ok<AGROW>

% Optional: de-duplicate (in case you re-run same pair); key = (pair,psi,N) etc.
[~, ia] = unique([ ...
    string(masterTable.pair_Jcap1) + "|" + ...
    string(masterTable.pair_Jcap2) + "|" + ...
    string(masterTable.psi) + "|" + ...
    string(masterTable.N)], 'stable');
masterTable = masterTable(ia,:);

save(matMaster, 'masterTable', '-v7');


fprintf('Saved sensitivity exports:\n  %s\n  %s\n(also appended %s)\n', ...
    csvPair, matPair, csvMaster);

% --- quick visuals --------------------------------------------------------
% (a) Z_norm vs φ colored by ψ
hZphi = figure('Visible','on'); hold on; grid on;
set(findall(hZphi,'Type','axes'),'FontSize',28,'LineWidth',1.2);
cmap = lines(numel(psiU));
for kpsi = 1:numel(psiU)
    m = (idPsi == kpsi);
    scatter(phiABn(m), ZABn(m), 30, cmap(kpsi,:), 'filled', ...
        'DisplayName', sprintf('\\psi = %.3g', psiU(kpsi)));
end
xlabel('\phi'); ylabel('F_{min}/(k_{eff} J_{intrin}^2)');
title('Z_{norm} vs \phi (colored by \psi)');
legend('Location','best');
set(legend(),'FontSize',20);
savefig(hZphi, fullfile(pairFolder, ...
    sprintf('Znorm_vs_phi_colored_by_psi_%0.4f_%0.4f.fig', aVal, bVal)));
close(hZphi);

% (b) Collapse check (reuse Zcollapse we computed already)
hCol = figure('Visible','on'); hold on; grid on;
set(findall(hCol,'Type','axes'),'FontSize',28,'LineWidth',1.2);
for kpsi = 1:numel(psiU)
    m = (idPsi == kpsi);
    plot(phiABn(m), Zcollapse(m), '.', 'Color', cmap(kpsi,:), ...
        'DisplayName', sprintf('\\psi = %.3g', psiU(kpsi)));
end
xlabel('\phi'); ylabel('Z_{norm} / mean_{band} (collapse check)');
title('Collapsed Z_{norm} (per-\psi mean removed)');
legend('Location','best');
set(legend(),'FontSize',20);
savefig(hCol, fullfile(pairFolder, ...
    sprintf('Znorm_collapse_by_psi_%0.4f_%0.4f.fig', aVal, bVal)));
close(hCol);

end
% ==========================================================================



fprintf('Finished.  Results appended to:\n   %s\n',matFile);
end

function plot_fit_results_from_MAT(matFilePath)
% PLOT_FIT_RESULTS_FROM_MAT
%   Reads the MAT-file produced by post-processing (variable = fitTable) and
%   generates three figures:
%        1) a   vs mismatch Δ  (Δ = |Jcap2 – Jcap1|)
%        2) b   vs mismatch Δ
%        3) R²  vs mismatch Δ
%
%   Each point is annotated by "[low, high]" where  low = min(Jcap1,Jcap2)
%   and high = max(Jcap1,Jcap2).
%
%   Rows whose b == 1 are ignored (as in the TXT version).
%
% Usage:
%   plot_fit_results_from_MAT('…/postproc Fmin fit results.mat')

%% --------------------------------------------------------------------- 
% 1)  Load table
%% ---------------------------------------------------------------------
S = load(matFilePath);
if ~isfield(S,'fitTable')
    error('MAT file %s does not contain variable ''fitTable''.',matFilePath);
end
T = S.fitTable;

% Guard against empty / header-only table
if isempty(T)
    error('fitTable is empty.');
end

%% ---------------------------------------------------------------------
% 2)  Build numeric arrays + labels
%% ---------------------------------------------------------------------
delta   = abs(T.Jcap2 - T.Jcap1);                % mismatch Δ
aVals   = T.a;
bVals   = T.b;
r2Vals  = T.Rsq;
aErrVals = T.a_err;

% Skip rows where b == 1  (per original rule)
keep    = (bVals ~= 1);
delta   = delta(keep);
aVals   = aVals(keep);
bVals   = bVals(keep);
r2Vals  = r2Vals(keep);
lowVals = min(T.Jcap1(keep), T.Jcap2(keep));
highVals= max(T.Jcap1(keep), T.Jcap2(keep));
labels  = arrayfun(@(lo,hi)sprintf('[%.4g, %.4g]',lo,hi), ...
                    lowVals, highVals, 'UniformOutput',false);

% Sort by Δ
[deltaSorted, idx] = sort(delta);
aSorted  = aVals(idx);
bSorted  = bVals(idx);
r2Sorted = r2Vals(idx);
aErrSorted = aErrVals(idx);
labelSorted = labels(idx);

%% helper for annotation
addLabels = @(x,y,lbl) arrayfun( ...
    @(i) text(x(i), y(i), lbl{i}, ...
              'VerticalAlignment','bottom','HorizontalAlignment','right', ...
              'FontSize',8,'Color','k'), ...
    1:numel(x));

%% ---------------------------------------------------------------------
% 3)  Figure:   a  vs Δ
%% ---------------------------------------------------------------------
[folderPath,~,~] = fileparts(matFilePath);

h1 = figure('Visible','on');
set(h1,'DefaultAxesFontSize',28,'DefaultAxesLineWidth',1.2);
plot(deltaSorted, aSorted, 'rs-','LineWidth',2,'MarkerSize',8);
xlabel('J_{prot} mismatch \Delta (High - Low)', ...
       'FontSize',30,'Color','red','FontWeight','bold');
ylabel('\alpha', 'FontSize',40,'Color','red','FontWeight','bold');
title('\alpha Values vs. Mismatch \Delta (with \kappa)', ...
      'FontSize',30,'Color','red','FontWeight','bold');
grid on;
ax = gca;  ax.FontSize = 28;  ax.LineWidth = 1.5;
ylim([1 2]);     % <<< ADD THIS LINE
addLabels(deltaSorted, aSorted, labelSorted);
mean_alpha    = mean(aSorted,'omitnan');
mean_alphaErr = std(aSorted,'omitnan');

txt = {sprintf('$\\alpha = %.3f \\pm %.3f$', ...
               mean_alpha, mean_alphaErr)};

annotation('textbox',[0.15 0.75 0.3 0.12], ...
    'String',txt, ...
    'Interpreter','latex', ...
    'FontSize',35, ...
    'Color','red', ...
    'EdgeColor','none', ...
    'FontWeight','bold');

savefig(h1, fullfile(folderPath,'aValues_vs_Delta.fig'));
close(h1);

%% ---------------------------------------------------------------------
% 4)  Figure:   b  vs Δ
%% ---------------------------------------------------------------------
h2 = figure('Visible','on');
set(h2,'DefaultAxesFontSize',28,'DefaultAxesLineWidth',1.2);
plot(deltaSorted, bSorted, 'bo-','LineWidth',2,'MarkerSize',8);
xlabel('J_{prot} mismatch \Delta (High - Low))', ...
       'FontSize',30,'Color','red','FontWeight','bold');
ylabel('\beta', 'FontSize',40,'Color','red','FontWeight','bold');
title('\beta Values vs. Mismatch \Delta (with \kappa)', ...
      'FontSize',30,'Color','red','FontWeight','bold');
grid on;
ax = gca;  ax.FontSize = 28;  ax.LineWidth = 1.5;
addLabels(deltaSorted, bSorted, labelSorted);
savefig(h2, fullfile(folderPath,'bValues_vs_Delta.fig'));
close(h2);

% keep only valid/positive
x = deltaSorted(:);  y = bSorted(:);
mask = isfinite(x) & isfinite(y) & x>0 & y>0;
x = x(mask);  y = y(mask);

% logs
X = log(x);                 % any base is fine
Y = log(y);

% ===== REPLACE FROM:  % ---- slope-only fit on centered logs ...  =====
% ---- full log-log linear fit WITH intercept:  log(y) = a + b*log(x) ----
X = log(x);
Y = log(y);

% Fit Y = a + b*X
[fo, gof] = fit(X, Y, 'poly1');   % fo.p1 = b (slope = gamma), fo.p2 = a (intercept = logC)
gamma = fo.p1;
logC  = fo.p2;
C     = exp(logC);
R2    = gof.rsquare;

% Confidence intervals (optional, nice to report)
ci = confint(fo);
gamma_CI = ci(:,1);
logC_CI  = ci(:,2);

% Line to plot (log-space)
xfit = linspace(min(X), max(X), 200);
yfit = gamma * xfit + logC;

% --- plot in log-log (linear axes are X=log(x), Y=log(y)) ---
h2 = figure('Visible','on');
set(h2,'DefaultAxesFontSize',28,'DefaultAxesLineWidth',1.2);
plot(X, Y, 'bo-','LineWidth',2,'MarkerSize',8); hold on
plot(xfit, yfit, 'r-','LineWidth',2);
grid on;
ax = gca;  ax.FontSize = 28;  ax.LineWidth = 1.5;
xlabel('log(\Delta)','FontSize',30,'Color','red','FontWeight','bold');
ylabel('log(\beta)','FontSize',30,'Color','red','FontWeight','bold');
title(sprintf('log–log fit: log\\beta = logC + \\gamma log\\Delta  |  \\gamma=%.3f  C=%.3g  R^2=%.3f', ...
              gamma, C, R2), ...
      'FontSize',22,'Color','red','FontWeight','bold');
addLabels(X, Y, labelSorted(mask));
legend({'data','fit: Y = a + bX'}, 'Location','best');
set(legend(),'FontSize',20);

savefig(h2, fullfile(folderPath,'bValues_vs_Delta_LOGLOG_withIntercept.fig'));
close(h2);

% --- OPTIONAL: plot in original (x,y) units with the power-law curve ---
% --- OPTIONAL: plot in original (x,y) units with the power-law curve ---
% ===== Recompute gamma in ORIGINAL SPACE =====
ftP = fittype('C*(x.^gamma)', ...
    'independent','x', ...
    'coefficients',{'C','gamma'});

optsP = fitoptions(ftP);
optsP.Lower = [0 0];
optsP.Upper = [Inf Inf];
optsP.StartPoint = [mean(y), 1];   % [C , gamma]

[fp, ~] = fit(x, y, ftP, optsP);

C     = fp.C;
gamma = fp.gamma;

% --- gamma error (95% CI)
try
    CIg = confint(fp);   % [C_low gamma_low ; C_high gamma_high]
    gamma_err = abs(CIg(2,2) - CIg(1,2))/2;
catch
    gamma_err = NaN;
end

h3 = figure('Visible','on');
set(h3,'DefaultAxesFontSize',28,'DefaultAxesLineWidth',1.2);

% data
plot(x, y, 'bo-','LineWidth',2,'MarkerSize',8); hold on

% model curve
xx = linspace(min(x), max(x), 200);
yy = C * xx.^gamma;
plot(xx, yy, 'r-','LineWidth',2);

% --- R^2 of the power fit in ORIGINAL (x,y) space ---
yhat = C * x.^gamma;
SSEp = sum((y - yhat).^2);
SSTp = sum((y - mean(y)).^2);
R2_power = 1 - SSEp/SSTp;
RMSE_power = sqrt(mean((y - yhat).^2));   % optional metric

grid on;
ax = gca;  ax.FontSize = 28;  ax.LineWidth = 1.5;
xlabel('J_{prot} mismatch \Delta (High - Low)', ...
       'FontSize',30,'Color','red','FontWeight','bold');
ylabel('\beta','FontSize',30,'Color','red','FontWeight','bold');

title('\beta Values vs. Mismatch \Delta', ...
      'FontSize',30,'Color','red','FontWeight','bold');
subtitle(sprintf('\\beta = C \\cdot \\Delta^{\\gamma}  |  \\gamma = %.3f \\pm %.3f   C = %.3g   (R^2 = %.3f, orig space)', ...
              gamma, gamma_err, C, R2_power), ...
      'FontSize',24,'Color','red','FontWeight','bold');


legend({'data','power-law fit'}, 'Location','best');
set(legend(),'FontSize',20);

% ---- add DeltaJvalues as labels on the points ----
addLabels(x, y, labelSorted(mask));
% ---- gamma textbox ----
txtG = {sprintf('$\\gamma = %.3f \\pm %.3f$', gamma, gamma_err)};

annotation('textbox',[0.15 0.75 0.3 0.12], ...
    'String',txtG, ...
    'Interpreter','latex', ...
    'FontSize',40, ...
    'Color', 'red', ...
    'EdgeColor','none', ...
    'FontWeight','bold', ...
    'Color','red');

savefig(h3, fullfile(folderPath,'bValues_vs_Delta_powerlaw_withIntercept.fig'));
close(h3);

% ===== REPLACE TO:  end of full log-log with intercept block =====


%% ---------------------------------------------------------------------
% 5)  Figure:   R² vs Δ
%% ---------------------------------------------------------------------
h3 = figure('Visible','on');
set(h3,'DefaultAxesFontSize',28,'DefaultAxesLineWidth',1.2);
plot(deltaSorted, r2Sorted, 'cd-','LineWidth',2,'MarkerSize',8);
xlabel('Mismatch \Delta (High - Low)', ...
       'FontSize',30,'Color','red','FontWeight','bold');
ylabel('R^2', 'FontSize',30,'Color','red','FontWeight','bold');
title('R^2 vs. Mismatch \Delta (with \kappa)', ...
      'FontSize',30,'Color','red','FontWeight','bold');
subtitle('Fmin - spherical model');
grid on;
ax = gca;  ax.FontSize = 28;  ax.LineWidth = 1.5;
addLabels(deltaSorted, r2Sorted, labelSorted);
savefig(h3, fullfile(folderPath,'R2Values_vs_Delta.fig'));
close(h3);

fprintf('Figures saved to %s\n', folderPath);
end


function [bestFit, bestGOF, psiClean, phiClean, ZClean, ...
          idxIn,   useLogFit, ZAllMaxBeforeLog] = ...
    fitWithMultipleStartPoints(aVal, bVal, ...
                               psiAB,  phiAB,  ZAB, ...
                               psiAA,  phiAA,  ZAA, ...
                               psiBB,  phiBB,  ZBB)
% =========================================================================
% Minimal change:
%   • If AA+BB are supplied (last six args) they are *only* used in the
%     fitting calls; every other step (pre-filter, residual trimming,
%     index bookkeeping) remains AB-only, so idxIn is still “within AB”.
%   • With the original 5-argument call, behaviour is 100 % identical.
% =========================================================================

% ---------- 0) fall-back when AA/BB not supplied --------------------------
if nargin < 11
    psiAA = [];  phiAA = [];  ZAA = [];
    psiBB = [];  phiBB = [];  ZBB = [];
end

% ---------- AB bookkeeping stays exactly the same -------------------------
origIdx = (1:numel(psiAB)).';        % AB-only indices

% ---------- assemble AB arrays for the *existing* workflow ---------------
psiAll = psiAB;   phiAll = phiAB;   ZAll = ZAB;   % everything below
                                                  % still talks “AB”

% ---------- 1) quick pre-filter on AB (unchanged) -------------------------
pctCut      = 99.9;
zCutValue   = prctile(ZAll, pctCut);
keep        = ZAll <= zCutValue;

psiAll  = psiAll(keep);
phiAll  = phiAll(keep);
ZAll    = ZAll  (keep);
origIdx = origIdx(keep);

% ---------- create combined arrays *only* for the fit --------------------
psiFit = [psiAll ; psiAA ; psiBB];
phiFit = [phiAll ; phiAA ; phiBB];
ZFit   = [ZAll   ; ZAA   ; ZBB  ];

% -- NEW: clamp Z to [0,1] -------------------------------------------
Zscale = max(ZFit);                % remember original scale
ZAll  = ZAll  / Zscale;
ZAB  = ZAB  / Zscale; 
ZAA   = ZAA   / Zscale;
ZBB   = ZBB   / Zscale;
ZFit  = ZFit  / Zscale;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>  raw-data scatter (AB + optional AA,BB)  – auto-closes 0.5 s  >>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
hFig = figure('Visible','on');  hold on;  grid on;  view(3);

% AB
scatter3(psiAB,phiAB,ZAB,45,'k','d','filled',...
         'DisplayName','AB raw');

% AA
if ~isempty(psiAA)
    scatter3(psiAA,phiAA,ZAA,40,'r','o','filled',...
             'DisplayName','AA raw');
end

% BB
if ~isempty(psiBB)
    scatter3(psiBB,phiBB,ZBB,40,'b','s','filled',...
             'DisplayName','BB raw');
end

xlabel('\psi'); ylabel('\phi'); zlabel('Z');
title('Raw data passed to the fit');
legend('Location','best');
set(legend(),'FontSize',20);
drawnow;
pause(0.5);
close(hFig);
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% ------------- log trick & guesses (unchanged, uses AB) ------------------
useLogFit = false;
ZAll(ZAll<=0)     = eps;
phiAll(phiAll<=0) = eps;
ZAllMaxBeforeLog = Zscale;

logPhi   = log(phiAll);
logZ     = log(ZAll);
a_guess  = max(0, polyfit(logPhi, logZ, 1));   % slope

b_guess1 = 1;
b_guess2 = mean(ZAll);

startPoints = [ ...
    a_guess(1), b_guess1;
    a_guess(1), b_guess2;
    0.5*a_guess(1), b_guess1;
    1.5*a_guess(1), b_guess2;
    2.0*a_guess(1), 0.5*b_guess1;
    0.8*a_guess(1), 2.0*b_guess2 ];

% ------------- model & options (unchanged) --------------------------------

if aVal < bVal
    deltaJ = bVal - aVal;
else
    deltaJ = aVal - bVal;
end

Kmem = 1;
ft = fittype(sprintf('(%g/2)*((psi.*(1-psi))).*(phi.^a).*b',Kmem), ...
             'independent',{'psi','phi'}, ...
             'coefficients',{'a','b'});
opts = fitoptions(ft);
opts.Robust      = 'Bisquare';
opts.Lower       = [0 0];
opts.Upper       = [5 Inf];
opts.MaxFunEvals = 1e4;
opts.MaxIter     = 1e4;

% ------------- 2) FIRST FIT (now uses AB+AA+BB) ---------------------------
ZFitTrans = ZFit;
bestR2 = -Inf;  bestFit = [];  bestGOF = [];
for k = 1:size(startPoints,1)
    opts.StartPoint = startPoints(k,:);
    try
        [f,g] = fit([psiFit,phiFit], ZFitTrans, ft, opts);
        if g.rsquare > bestR2
            bestR2 = g.rsquare;  bestFit = f;  bestGOF = g;
        end
    end
end
if isempty(bestFit), error('All fits failed in the initial pass.'); end

% ------------- rescale b if log-trick (unchanged) -------------------------
c = coeffvalues(bestFit);
aFit = c(1);  bFit = c(2);


% ------------- 3) residuals & optional trimming (AB-only, unchanged) ------
ZModel = (Kmem/2) * bFit * ...                        % scalar factors
         ((psiAll .* (1-psiAll))) ...              % (ψ·(1-ψ))²
         .* (phiAll .^ aFit);                      % φ^a;                             % divide by (1-φ)
resid     = abs(ZAll - ZModel);

% the rest of the original code (trimming, etc.) follows ***unchanged***
% --------------------------------------------------------------------------
rRange  = max(resid) - min(resid);
rThresh = 1e-5;
if rRange >= rThresh
    maxIters     = 1;
    trimFraction = 1e-6;
    minPoints    = 10;

    for iter = 1:maxIters
        % iterate *exactly as before* on AB data only
        ZTransAB = ZAll;
        bestR2 = -Inf; bestFit = []; bestGOF = [];
        for k = 1:size(startPoints,1)
            opts.StartPoint = startPoints(k,:);
            try
                [f,g] = fit([psiAll,phiAll], ZTransAB, ft, opts);
                if g.rsquare > bestR2
                    bestR2= g.rsquare; bestFit = f; bestGOF = g;
                end
            end
        end
        c = coeffvalues(bestFit);
        % if useLogFit
        %     bestFit.b = bestFit.b / 1e6;      % <<< ADD THIS LINE
        %     aFit = c(1); bFit = c(2)/1e6;
        % else
            aFit = c(1); bFit = c(2);
        % end
        ZModel = (Kmem/2) * bFit * ...                        % scalar factors
                 ((psiAll .* (1-psiAll))) ...              % (ψ·(1-ψ))²
                 .* (phiAll .^ aFit);                      % φ^a;                             % divide by (1-φ)
        resid  = abs(ZAll - ZModel);
        N = numel(ZAll);
        keepN = max(minPoints, round((1-trimFraction)*N));
        [~,ord] = sort(resid,'ascend');
        keep    = ord(1:keepN);
        psiAll  = psiAll(keep);
        phiAll  = phiAll(keep);
        ZAll    = ZAll  (keep);
        origIdx = origIdx(keep);
    end
end


% ------------- 4) outputs -------------------------------------------------
bestFit.b = bestFit.b * Zscale;
psiClean = psiAll;
phiClean = phiAll;
ZClean   = ZAll * Zscale;
idxIn    = origIdx;               % still AB indices only
end



% ========================================================================
function addSymmetricScatter(d, labelFmt, JcapVal, protSet, cMap, mkr)
% Plot coloured AA or BB points, one legend entry per protein colour.
if isempty(d.phi), return, end
for p = protSet(:).'
    idx = d.prot == p;
    if ~any(idx), continue, end
    rgbC = cMap(p,:);
    lbl  = sprintf('%s  proteins = %d', sprintf(labelFmt,JcapVal,JcapVal), p);
    scatter3(d.phi(idx), d.psi(idx), d.Z(idx), ...
             40, rgbC, mkr, 'filled', 'DisplayName', lbl);
end
end

function addABInlierScatter(phiIn, psiIn, ZIn, protIn, labelFmt, aVal, bVal, protSet, cMap, orientationFlag)
% orientationFlag: 1 = AB  → use AA marker ('^')
%                  2 = BA  → use BB marker ('s')

    % pick marker shape according to orientation
    switch orientationFlag
        case 1
            markerShape = '^';   % same as AA
        case 2
            markerShape = 's';   % same as BB
        otherwise
            markerShape = 'o';   % fallback, should not happen
    end

    for p = protSet(:).'
        m = (protIn == p);
        if ~any(m), continue; end

        scatter3(phiIn(m), psiIn(m), ZIn(m), ...
                 45, cMap(p,:), markerShape, 'filled', ...
                 'DisplayName', sprintf([labelFmt '  proteins = %d'], aVal, bVal, p));
    end
end



function [phiAB_filt, psiAB_filt, Z_AB_filt, protAB_filt, midJcapAB_filt, cordAB_filt] = filter_AB_by_PSI_fit_quality(phiAB, psiAB, Z_AB, protAB, ...
                                  tolPsi, minRSquare, minPointsPerPsi)

% FILTER_AB_BY_PSI_FIT_QUALITY
% -------------------------------------------------------------------------
% Filters AB data by evaluating the fit quality (R²) of Z vs φ at each ψ slice.
% Discards entire ψ bands where quadratic fit R² is below minRSquare.
%
% Inputs:
%   phiAB, psiAB, Z_AB, protAB, midJcapAB, cordAB  : AB data arrays
%   tolPsi        : tolerance for grouping ψ (default: 1e-6)
%   minRSquare    : minimum R² to keep a ψ slice (default: 0.95)
%   minPointsPerPsi : minimum points per ψ slice to attempt fitting (default: 5)
%
% Outputs:
%   Same inputs but filtered
% -------------------------------------------------------------------------

if nargin < 7, tolPsi = 1e-6; end
if nargin < 8, minRSquare = 0.95; end
if nargin < 9, minPointsPerPsi = 5; end

[uniquePsi,~,idPsi] = uniquetol(psiAB, tolPsi);

keepMask = false(size(psiAB));   % initialize as all false

for k = 1:numel(uniquePsi)
    idxThisPsi = idPsi == k;
    
    if sum(idxThisPsi) < minPointsPerPsi
        continue
    end
    
    % Fit quadratic
    phi_k = phiAB(idxThisPsi);
    Z_k   = Z_AB(idxThisPsi);
    
    p = polyfit(phi_k, Z_k, 2);
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
phiAB_filt     = phiAB(keepMask);
psiAB_filt     = psiAB(keepMask);
Z_AB_filt      = Z_AB(keepMask);
protAB_filt    = protAB(keepMask);
midJcapAB_filt = midJcapAB(keepMask);
cordAB_filt    = cordAB(keepMask);

end


function [phi_out, Z_out, prot_out] = filter_by_protein_fit_quality(phi, Z, prot, minRSquare, minPointsPerProt)

    if nargin < 4, minRSquare = 0.95; end
    if nargin < 5, minPointsPerProt = 5; end

    phi_out = [];
    Z_out = [];
    prot_out = [];

    uniqueProt = unique(prot);

    for i = 1:numel(uniqueProt)
        p = uniqueProt(i);
        idx = (prot == p);

        % Extract and clean this protein's data
        phi_p = phi(idx);
        Z_p   = Z(idx);

        validMask = isfinite(phi_p) & isfinite(Z_p);
        phi_p = phi_p(validMask);
        Z_p   = Z_p(validMask);

        if numel(Z_p) < minPointsPerProt
            fprintf('prot %d: skipped (N = %d < %d)\n', p, numel(Z_p), minPointsPerProt);
            continue;
        end

        % Check for zero variance
        if std(phi_p) == 0 || std(Z_p) == 0
            R2 = NaN;
        else
            % Linear fit
            pfit = polyfit(phi_p, Z_p, 2);
            Z_fit = polyval(pfit, phi_p);

            % Compute R² safely
            SS_res = sum((Z_p - Z_fit).^2);
            SS_tot = sum((Z_p - mean(Z_p)).^2);
            if SS_tot == 0
                R2 = NaN;
            else
                R2 = 1 - SS_res / SS_tot;
            end
        end

        fprintf('prot %d: N = %d → R² = %.4g\n', p, numel(Z_p), R2);

        % Keep only if good
        if R2 >= minRSquare
            phi_out  = [phi_out ; phi_p];
            Z_out    = [Z_out   ; Z_p];
            prot_out = [prot_out; repmat(p, numel(phi_p), 1)];
        end
    end
end



function r = tiedrank_fallback(x)
% Tied-rank without Statistics Toolbox (average ranks for ties).
[xs, ord] = sort(x(:));
n = numel(xs);
ranks = zeros(n,1);
i = 1;
while i <= n
    j = i;
    while j < n && xs(j+1) == xs(i)
        j = j + 1;
    end
    ranks(i:j) = (i + j) / 2;   % average rank for the tie block
    i = j + 1;
end
% put back to original order
r = zeros(n,1);
r(ord) = ranks;
end

%% *************************************************************************
function build_deltaJcap_scans(output_path)
% BUILD_DELTAJCAP_SCANS
% Aggregates AB-inlier caches from pair_* folders under output_path and,
% for each ψ, creates a figure with:
%   Z = Fmin,  X = phi,  Y = ΔJcap
% Then fits the custom surface:
%   Fmin = (Kmem/2)*(ψ(1-ψ))^2 * (phi^α) * (ΔJcap)^β
% Saves .fig + points CSV + per-ψ fit .mat (+ a running fit_summary.csv).
%
% Requires that each pair folder contains:
%   pair_AB_inliers_*.mat  with struct 'pairData' having fields:
%   Jcap1, Jcap2, deltaJcap, phi, psi, Fmin
%
% Folder layout produced:
%   <output_path>/deltaJcaps/PSI <ψ_key>/
%       Z_vs_phi_vs_deltaJcap_PSI_<ψ_key>.fig
%       points_PSI_<ψ_key>.csv
%       fit_PSI_<ψ_key>.mat
%   and a running:
%       <output_path>/deltaJcaps/fit_summary.csv

% --- sanity
if ~exist(output_path,'dir')
    error('Output path does not exist: %s', output_path);
end

pairDirs = dir(fullfile(output_path,'pair_*'));
if isempty(pairDirs)
    warning('No pair_* folders found under %s', output_path);
    return;
end

% ---- collect all AB-inlier points across pairs ---------------------------
PHI = []; PSI = []; FMIN = []; DELTA = [];
for d = 1:numel(pairDirs)
    if ~pairDirs(d).isdir, continue; end
    pairDir = fullfile(pairDirs(d).folder, pairDirs(d).name);
    M = dir(fullfile(pairDir, 'pair_AB_inliers_*.mat'));
    if isempty(M), continue; end
    for m = 1:numel(M)
        S = load(fullfile(M(m).folder, M(m).name));
        if ~isfield(S,'pairData'), continue; end
        pd = S.pairData;

        good = isfinite(pd.phi) & isfinite(pd.psi) & isfinite(pd.Fmin);
        if any(good)
            PHI   = [PHI  ; pd.phi(good)];
            PSI   = [PSI  ; pd.psi(good)];
            FMIN  = [FMIN ; pd.Fmin(good)];
            DELTA = [DELTA; repmat(pd.deltaJcap, nnz(good), 1)];
        end
    end
end

if isempty(PHI)
    warning('No AB inlier data found to build ΔJcap figures.');
    return;
end

% ---- collision-proof ψ grouping & naming (fixed decimals) ----------------
NDIG = 6;                                   % folder/key precision for ψ
psiRounded = round(PSI, NDIG);
psiKeys    = arrayfun(@(v) sprintf(['%0.' num2str(NDIG) 'f'], v), psiRounded, 'uni', false);
uKeys      = unique(psiKeys, 'stable');

deltaRoot = fullfile(output_path, 'deltaJcaps');
if ~exist(deltaRoot,'dir'), mkdir(deltaRoot); end

% Prepare global fit summary CSV (append; write header once)
summaryCSV  = fullfile(deltaRoot, 'fit_summary.csv');
writeHeader = ~isfile(summaryCSV);

for k = 1:numel(uKeys)
    key  = uKeys{k};                           % e.g., '0.833333'
    mask = strcmp(psiKeys, key);
    if ~any(mask), continue; end

    psiFolder = fullfile(deltaRoot, ['PSI ' key]);  % "PSI 0.833333"
    if ~exist(psiFolder,'dir'), mkdir(psiFolder); end

    x_phi   = PHI(mask);
    y_delta = DELTA(mask);
    z_fmin  = FMIN(mask);

    % ---------- custom surface fit: Fmin = Cpsi * (phi^a) * (delta^b)
    % domain guards
    maskOK = isfinite(x_phi) & isfinite(y_delta) & isfinite(z_fmin) ...
           & (x_phi > 0) & (x_phi < 1) & (y_delta > 0) & (z_fmin > 0);
    x_phi   = x_phi(maskOK);
    y_delta = y_delta(maskOK);
    z_fmin  = z_fmin(maskOK);
    if isempty(x_phi), continue; end

    % numeric ψ from the key; Kmem default 1 (change here if needed)
    psiNum = str2double(key);
    Kmem   = 1;
    Cpsi   = (Kmem/2) * (psiNum*(1-psiNum))^2;

    % fit type and options
    ft = fittype( sprintf('%g * ((phi.^a)) .* (delta.^b)', Cpsi), ...
                  'independent', {'phi','delta'}, 'coefficients', {'a','b'} );
    opts = fitoptions(ft);
    opts.Robust      = 'Bisquare';
    opts.Lower       = [0 0];
    opts.Upper       = [5 Inf];
    opts.MaxFunEvals = 1e4;
    opts.MaxIter     = 1e4;

    % multistart
    aStarts = [0.2 0.6 1.0 1.5 2.0 2.5 3.0];
    bStarts = [0.2 0.6 1.0 1.5 2.0 2.5 3.0];
    bestGOF = struct('rsquare', -Inf);
    bestFit = [];

    if ~isempty(x_phi)
        for aa = aStarts
            for bb = bStarts
                opts.StartPoint = [aa, bb];
                try
                    [f, g] = fit([x_phi, y_delta], z_fmin, ft, opts);
                    if g.rsquare > bestGOF.rsquare
                        bestGOF = g;  bestFit = f;
                    end
                catch
                    % try next start
                end
            end
        end
    end

    % ---- plot scatter + fitted surface (if fit succeeded)
    h = figure('Visible','on'); hold on; grid on; view(3);
    scatter3(x_phi, y_delta, z_fmin, 30, 'filled');

    eqStr = ['Fit Equation: F_{min} = (K_{mem}/2) (\psi(1-\psi))^{2} \cdot ' ...
             '(\phi^{\alpha}) \cdot (\Delta J_{cap})^{\beta}'];

    if ~isempty(bestFit)
        phiMesh   = linspace(min(x_phi),   max(x_phi),   40);
        deltaMesh = linspace(min(y_delta), max(y_delta), 40);
        [PhiM, DelM] = meshgrid(phiMesh, deltaMesh);
        Zhat = feval(bestFit, PhiM, DelM);
        surf(PhiM, DelM, Zhat, 'EdgeColor','none', 'FaceAlpha',0.35);

        c = coeffvalues(bestFit);           % [a b]
        aHat = c(1);  bHat = c(2);
        fitInfoStr = sprintf('\\alpha = %.3g, \\beta = %.3g, R^2 = %.3f', aHat, bHat, bestGOF.rsquare);
        titleStr = sprintf('Z=F_{min}, X=\\phi, Y=\\Delta J_{cap} (\\psi = %s) | %s', key, fitInfoStr);
    else
        titleStr = sprintf('Z=F_{min}, X=\\phi, Y=\\Delta J_{cap} (\\psi = %s)', key);
    end

    xlabel('\phi', 'FontSize',14, 'FontWeight','bold');
    ylabel('\Delta J_{cap} = |J_2 - J_1|', 'FontSize',14, 'FontWeight','bold');
    zlabel('F_{min}', 'FontSize',14, 'FontWeight','bold');
    title(titleStr, 'FontSize',14, 'FontWeight','bold');
    subtitle({eqStr}, 'Interpreter','tex', 'FontSize',20, ...
             'Color','black', 'FontWeight','bold');

    savefig(h, fullfile(psiFolder, sprintf('Z_vs_phi_vs_deltaJcap_PSI_%s.fig', key)));
    close(h);

    % ---- save points for this ψ
    T = table(x_phi, y_delta, z_fmin, 'VariableNames', {'phi','deltaJcap','Fmin'});
    writetable(T, fullfile(psiFolder, sprintf('points_PSI_%s.csv', key)));

    % ---- save per-ψ fit results + append to summary
    if ~isempty(bestFit)
        c = coeffvalues(bestFit);
        aHat = c(1);  bHat = c(2);
        FitResult = struct();
        FitResult.psi       = psiNum;
        FitResult.Kmem      = Kmem;
        FitResult.Cpsi      = Cpsi;
        FitResult.alpha     = aHat;
        FitResult.beta      = bHat;
        FitResult.Rsq       = bestGOF.rsquare;
        FitResult.modelTex  = 'F_{min}=(K_{mem}/2)(\psi(1-\psi))\cdot(\phi^{\alpha})\cdot(\Delta J_{cap})^{\beta}';
        save(fullfile(psiFolder, sprintf('fit_PSI_%s.mat', key)), 'FitResult');

        % summary CSV (header once)
        fid = fopen(summaryCSV, 'a');
        if fid > 0
            if writeHeader
                fprintf(fid, 'psi,alpha,beta,Rsq,Kmem,Cpsi,points\n');
                writeHeader = false;
            end
            fprintf(fid, '%.10g,%.6g,%.6g,%.6g,%.6g,%.6g,%d\n', ...
                psiNum, aHat, bHat, bestGOF.rsquare, Kmem, Cpsi, numel(x_phi));
            fclose(fid);
        end
    end
end

fprintf('ΔJcap figures + fits written to:\n  %s\n', deltaRoot);
end




%% *****************************************************************************
function plot_deltaJcap_fit_results(output_path)
% PLOT_DELTAJCAP_FIT_RESULTS
% Reads <output_path>/deltaJcaps/fit_summary.csv (written by
% build_deltaJcap_scans) and generates three figures:
%   1) alpha vs psi
%   2) beta  vs psi
%   3) R^2   vs psi
%
% Usage:
%   plot_deltaJcap_fit_results(output_path)

%% ---------------------------------------------------------------------
% 1) Locate & read summary CSV
%% ---------------------------------------------------------------------
deltaRoot = fullfile(output_path, 'deltaJcaps');
csvPath   = fullfile(deltaRoot, 'fit_summary.csv');

if ~isfile(csvPath)
    error('File not found: %s (run build_deltaJcap_scans first)', csvPath);
end

T = readtable(csvPath);
req = {'psi','alpha','beta','Rsq'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames)
        error('Column "%s" not found in %s', req{i}, csvPath);
    end
end

% guard against empties / NaNs
mask = isfinite(T.psi) & isfinite(T.alpha) & isfinite(T.beta) & isfinite(T.Rsq);
T = T(mask, :);
if isempty(T)
    error('No valid rows in %s', csvPath);
end

% sort by psi (ascending)
[Tpsi, idx] = sort(T.psi);
Ta = T.alpha(idx);
Tb = T.beta(idx);
Tr = T.Rsq(idx);

% labels (formatted ψ) for annotations
labels = arrayfun(@(v) sprintf('\\psi = %.6g', v), Tpsi, 'uni', false);

% helper to annotate points
addLabels = @(x,y,lbl) arrayfun( ...
    @(i) text(x(i), y(i), lbl{i}, ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','right', ...
        'FontSize',8, 'Color','k'), ...
    1:numel(x));

%% ---------------------------------------------------------------------
% 2) α vs ψ
%% ---------------------------------------------------------------------
h1 = figure('Visible','on');
plot(Tpsi, Ta, 'rs-','LineWidth',2,'MarkerSize',8);
grid on;
xlabel('\psi',   'FontSize',30,'Color','red','FontWeight','bold');
ylabel('\alpha', 'FontSize',30,'Color','red','FontWeight','bold');
title('\alpha vs. \psi (from F_{min} fit: \Delta J_{cap}^{\beta})', ...
      'FontSize',26,'Color','red','FontWeight','bold');
addLabels(Tpsi, Ta, labels);
savefig(h1, fullfile(deltaRoot, 'alpha_vs_psi.fig'));
close(h1);

%% ---------------------------------------------------------------------
% 3) β vs ψ
%% ---------------------------------------------------------------------
h2 = figure('Visible','on');
plot(Tpsi, Tb, 'bo-','LineWidth',2,'MarkerSize',8);
grid on;
xlabel('\psi',  'FontSize',30,'Color','red','FontWeight','bold');
ylabel('\beta', 'FontSize',30,'Color','red','FontWeight','bold');
title('\beta vs. \psi (exponent on \Delta J_{cap})', ...
      'FontSize',26,'Color','red','FontWeight','bold');
addLabels(Tpsi, Tb, labels);
savefig(h2, fullfile(deltaRoot, 'beta_vs_psi.fig'));
close(h2);

%% ---------------------------------------------------------------------
% 4) R^2 vs ψ
%% ---------------------------------------------------------------------
h3 = figure('Visible','on');
plot(Tpsi, Tr, 'cd-','LineWidth',2,'MarkerSize',8);
grid on;
xlabel('\psi', 'FontSize',30,'Color','red','FontWeight','bold');
ylabel('R^2',  'FontSize',30,'Color','red','FontWeight','bold');
title('R^2 vs. \psi (goodness of fit)', ...
      'FontSize',26,'Color','red','FontWeight','bold');
addLabels(Tpsi, Tr, labels);
savefig(h3, fullfile(deltaRoot, 'R2_vs_psi.fig'));
close(h3);

fprintf('Per-\\psi fit figures saved to %s\n', deltaRoot);
end
