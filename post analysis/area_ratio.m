close all
clear all
clc

% ====== paths ======
output_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon\V14\my_post_proc_figures\AreaRatio\';
postStructmat_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon\V14\my_post_proc_figures\postProcessingStruct_Intrinsic.mat';

% ====== Jcap grids to sweep ======
Jcap1 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
Jcap2 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];

% ====== run ======
run_for_all_Jcap_combinations_AR(Jcap1, Jcap2, postStructmat_path, output_path);
build_symmetricJcap_scan_AR(postStructmat_path, output_path);
% build_deltaJcap_scans_AR(output_path);

%% =======================================================================
function run_for_all_Jcap_combinations_AR(Jcap1, Jcap2, postStructmat_path, output_path)

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    if isempty(Jcap1) || isempty(Jcap2)
        error('Jcap1 and Jcap2 must not be empty');
    end

    processedPairs = containers.Map();  % keys like '0.0010_0.1500'

    for j = 1:length(Jcap1)
        for k = 1:length(Jcap2)
            N1 = Jcap1(j);
            N2 = Jcap2(k);

            if abs(N1 - N2) < 1e-12
                fprintf('Skipping identical pair: [%.4f, %.4f]\n', N1, N2);
                continue;
            end

            % unordered key to avoid mirrored duplicates
            key = sprintf('%.10f_%.10f', min(N1, N2), max(N1, N2));
            if isKey(processedPairs, key)
                fprintf('Skipping mirrored duplicate: [%.4f, %.4f]\n', N1, N2);
                continue;
            end
            processedPairs(key) = true;

            fprintf('Processing pair [%.4f, %.4f]\n', N1, N2);
            try
                postproc_area_ratio_selected_pair_struct(postStructmat_path, output_path, N1, N2);
            catch ME
                warning('Failed processing pair [%.4f, %.4f]: %s', N1, N2, ME.message);
            end
        end
    end
end

%% =======================================================================
function postproc_area_ratio_selected_pair_struct(structPath, output_path, Jcap1, Jcap2)
% POSTPROC_AREA_RATIO_SELECTED_PAIR_STRUCT
% - Loads the post-processing struct and collects data for the unordered pair {Jcap1,Jcap2}.
% - Pulls AB rows (a≠b), plus [a,a] and [b,b] if present.
% - Z-values are shown as **percent deviation from sphere**:
%       100 * (real - analytic) / analytic
% - **No filtering. No fitting.** Just plot 3-D scatters/surface in (phi, psi, % deviation).
% - Also saves all AB points for this pair for later ΔJcap aggregation
%   (keeps original area_ratio and adds area_ratio_pct).

tol  = 1e-6;
aVal = Jcap1;
bVal = Jcap2;
if abs(aVal-bVal) < tol
    error('Jcap1 and Jcap2 must be different.');
end

% -- load struct --
S = load(structPath);
fn = fieldnames(S);
if numel(fn) ~= 1 || ~isstruct(S.(fn{1}))
    error('File %s must contain ONE struct variable.', structPath);
end
S = S.(fn{1});

% detect area-ratio field
hasAR_Fmin = isfield(S, 'area_ratio_at_Fmin');
hasAR      = isfield(S, 'area_ratio');
if ~hasAR_Fmin && ~hasAR
    error('Struct has neither "area_ratio_at_Fmin" nor "area_ratio".');
end

% -- split AB / AA / BB (no filtering) --
phiAB = [];  psiAB = [];  Z_AB = [];
dataAA = struct('phi',[],'psi',[],'Z',[]);
dataBB = struct('phi',[],'psi',[],'Z',[]);

for k = 1:numel(S)
    jcaps = unique(S(k).Jcap(:));
    jcaps = sort(jcaps);

    % choose Z source for this entry
    if hasAR_Fmin
        Zi = S(k).area_ratio_at_Fmin(:);
    else
        Zi = S(k).area_ratio(:);
    end

    % AB (mixed)
    if numel(jcaps)==2 && ismembertol(aVal,jcaps,tol) && ismembertol(bVal,jcaps,tol)
        phiAB = [phiAB; S(k).phi(:)];
        psiAB = [psiAB; S(k).psi(:)];
        Z_AB  = [Z_AB ; Zi];
        continue
    end
    % AA
    if numel(jcaps)==1 && abs(jcaps - aVal) < tol
        dataAA.phi = [dataAA.phi; S(k).phi(:)];
        dataAA.Z   = [dataAA.Z  ; Zi];
        continue
    end
    % BB
    if numel(jcaps)==1 && abs(jcaps - bVal) < tol
        dataBB.phi = [dataBB.phi; S(k).phi(:)];
        dataBB.Z   = [dataBB.Z  ; Zi];
    end
end

if isempty(phiAB)
    error('Pair [%.4f, %.4f] not found in %s', aVal, bVal, structPath);
end

% hard-code ψ for symmetric sets so they sit at edges
if aVal < bVal
    dataAA.psi = ones(size(dataAA.phi));   % AA at ψ=1 when a is the smaller
    dataBB.psi = zeros(size(dataBB.phi));  % BB at ψ=0
else
    dataAA.psi = zeros(size(dataAA.phi));
    dataBB.psi = ones(size(dataBB.phi));
end

% --- convert area-ratio to **percent deviation** from sphere --------------
% (pure addition; original Z arrays kept for saving too)
Z_AB_pct     = 100 .* Z_AB;
dataAA.Z_pct = 100 .* dataAA.Z;
dataBB.Z_pct = 100 .* dataBB.Z;

% --- plotting (phi, psi, % deviation) ---
if ~exist(output_path,'dir'), mkdir(output_path); end
pairTag    = sprintf('%.4f_%.4f', aVal, bVal);   % order-preserving tag
pairFolder = fullfile(output_path, ['pair_' pairTag]);
if ~exist(pairFolder,'dir'), mkdir(pairFolder); end

hFig = figure('Visible','on'); hold on; grid on; view(3);

% === Half-transparent surface over (phi, psi) =============================
% Gather points to interpolate (AB + AA + BB anchor edges) in PERCENT
phiFit = [phiAB(:); dataAA.phi(:); dataBB.phi(:)];
psiFit = [psiAB(:); dataAA.psi(:); dataBB.psi(:)];
Zfit   = [Z_AB_pct(:); dataAA.Z_pct(:); dataBB.Z_pct(:)];

maskValid = isfinite(phiFit) & isfinite(psiFit) & isfinite(Zfit);
phiFit = phiFit(maskValid);  psiFit = psiFit(maskValid);  Zfit = Zfit(maskValid);

% Require at least 3 unique samples
[XYu, uniqIdx] = unique([phiFit psiFit], 'rows');
phiFit = XYu(:,1); psiFit = XYu(:,2); Zfit = Zfit(uniqIdx);
if numel(Zfit) >= 3
    % Interpolated surface (smooth) — NOT a fit
    nGrid = 120;
    phiL = min(phiFit); phiU = max(phiFit);
    psiL = 0;           psiU = 1;
    [PhiG, PsiG] = meshgrid(linspace(phiL,phiU,nGrid), linspace(psiL,psiU,nGrid));

    % Smooth interpolation that passes through data
    method = 'v4';  % biharmonic spline
    ZG = griddata(phiFit, psiFit, Zfit, PhiG, PsiG, method);

    % Keep only the data convex hull (outside becomes NaN)
    K  = convhull(phiFit, psiFit);
    in = inpolygon(PhiG, PsiG, phiFit(K), psiFit(K));
    ZG(~in) = NaN;

    % Draw half-transparent surface under the points
    set(gcf,'Renderer','opengl');  % transparency in PNGs
    sH = surf(PhiG, PsiG, ZG, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.50, ...
        'DisplayName','Interpolated surface');
    shading interp;
    material dull; lighting gouraud; camlight headlight;
end

% === Color scale bar for the surface (percent) ============================
if exist('sH','var') && isgraphics(sH)
    % Set color limits from the surface (ignoring NaNs)
    zSurfMin = min(ZG(:), [], 'omitnan');
    zSurfMax = max(ZG(:), [], 'omitnan');
    if isfinite(zSurfMin) && isfinite(zSurfMax) && zSurfMin < zSurfMax
        caxis([zSurfMin zSurfMax]);
    end

    % Add colorbar (scale bar) with percent tick labels
    cb = colorbar('Location','eastoutside');

    % ticks on LEFT side of the bar
    cb.YAxisLocation = 'left';

    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 22;

    % ONLY bottom + top ticks
    cb.Ticks = [zSurfMin, zSurfMax];

    % 2 decimals, percent formatting
    dec = 2;
    fmt = sprintf('%%.%df%%%%', dec);
    cb.TickLabels = {sprintf(fmt, zSurfMin), sprintf(fmt, zSurfMax)};

    set(gcf,'Renderer','opengl');  % transparency-friendly
end


% === Points (use percent for Z) ===========================================
% AB: all points (black diamonds)
scatter3(phiAB, psiAB, Z_AB_pct, 42, 'k', 'd', 'filled', ...
         'DisplayName', sprintf('[%.4f, %.4f]  AB', aVal, bVal));

% AA: red triangles (if present)
if ~isempty(dataAA.phi)
    scatter3(dataAA.phi, dataAA.psi, dataAA.Z_pct, ...
             50, 'r', '^', 'filled', ...
             'DisplayName', sprintf('[%.4f, %.4f]  AA', aVal, aVal));
end

% BB: blue squares (if present)
if ~isempty(dataBB.phi)
    scatter3(dataBB.phi, dataBB.psi, dataBB.Z_pct, ...
             50, 'b', 's', 'filled', ...
             'DisplayName', sprintf('[%.4f, %.4f]  BB', bVal, bVal));
end

annotation(gcf,'textbox',[0.02 0.35 0.10 0.30], ...   % [x y w h] (normalized)
    'String', {'$\frac{\left(\mathrm{micro\ area}-\mathrm{macro\ area}\right)}{\mathrm{macro\ area}}$', ...
'', ...
'$[\%]$'}, ...
    'Interpreter','latex', ...
    'FontSize',35, ...
    'FontWeight','bold', ...
    'Color','red', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');

xlabel('\phi', 'FontSize',45, 'Color','red', 'FontWeight','bold');
ylabel('\psi', 'FontSize',45, 'Color','red', 'FontWeight','bold');
set(gcf,'Renderer','opengl');

zlabel('Area deviation (%)', ...
       'Interpreter','tex', ...
       'FontSize',45, 'Color','red', 'FontWeight','bold');

title(sprintf(['Area deviation ($\\frac{\\mathrm{micro\\ area}-\\mathrm{macro\\ area}}{\\mathrm{macro\\ area}}\\;[\\%%]$) ' ...
               '- pair [%.4f, %.4f]'], aVal, bVal), ...
      'Interpreter','latex', ...
      'FontSize',50, 'Color','red', 'FontWeight','bold');

legend('Location','best');

savefig(hFig, fullfile(pairFolder, sprintf('AreaRatio_phi_psi_%s.fig', pairTag)));
saveas(hFig, fullfile(pairFolder, sprintf('AreaRatio_phi_psi_%s.png', pairTag)));
close(hFig);

% --- save AB points for ΔJcap aggregation (no filtering) ------------------
% Keep original ratio AND add the percent field
deltaJcap = abs(aVal - bVal);
pairData = struct( ...
    'Jcap1', aVal, ...
    'Jcap2', bVal, ...
    'deltaJcap',     deltaJcap, ...
    'phi',           phiAB(:), ...
    'psi',           psiAB(:), ...
    'area_ratio',    Z_AB(:), ...        % original unitless ratio
    'area_ratio_pct',Z_AB_pct(:) );      % percent deviation (added)

save(fullfile(pairFolder, sprintf('pair_AB_points_%s.mat', pairTag)), 'pairData', '-v7');

fprintf('Saved percent-deviation plot + points for pair [%g, %g]\n', aVal, bVal);
end


%% =======================================================================
function build_deltaJcap_scans_AR(output_path)
% BUILD_DELTAJCAP_SCANS_AR
% Aggregates AB point-caches from pair_* folders under output_path and,
% for each ψ (grouped by tolerance), creates a figure with:
%   Z = area_ratio,  X = φ,  Y = ΔJcap
% No fitting, no filtering — just scatter plots + CSVs.

% sanity
if ~exist(output_path,'dir')
    error('Output path does not exist: %s', output_path);
end

pairDirs = dir(fullfile(output_path,'pair_*'));
if isempty(pairDirs)
    warning('No pair_* folders found under %s', output_path);
    return;
end

% collect AB points across all pairs
PHI = []; PSI = []; AR = []; DELTA = [];
for d = 1:numel(pairDirs)
    if ~pairDirs(d).isdir, continue; end
    pairDir = fullfile(pairDirs(d).folder, pairDirs(d).name);
    M = dir(fullfile(pairDir, 'pair_AB_points_*.mat'));
    if isempty(M), continue; end
    for m = 1:numel(M)
        S = load(fullfile(M(m).folder, M(m).name));
        if ~isfield(S,'pairData'), continue; end
        pd = S.pairData;

        good = isfinite(pd.phi) & isfinite(pd.psi) & isfinite(pd.area_ratio);
        if any(good)
            PHI   = [PHI  ; pd.phi(good)];
            PSI   = [PSI  ; pd.psi(good)];
            AR    = [AR   ; pd.area_ratio(good)];
            DELTA = [DELTA; repmat(pd.deltaJcap, nnz(good), 1)];
        end
    end
end

if isempty(PHI)
    warning('No AB data found to build ΔJcap area-ratio figures.');
    return;
end

% group by ψ using fixed precision keys
NDIG = 6;
psiRounded = round(PSI, NDIG);
psiKeys    = arrayfun(@(v) sprintf(['%0.' num2str(NDIG) 'f'], v), psiRounded, 'uni', false);
uKeys      = unique(psiKeys, 'stable');

deltaRoot = fullfile(output_path, 'deltaJcaps_area_ratio');
if ~exist(deltaRoot,'dir'), mkdir(deltaRoot); end

for k = 1:numel(uKeys)
    key  = uKeys{k};                 % e.g., '0.833333'
    mask = strcmp(psiKeys, key);
    if ~any(mask), continue; end

    psiFolder = fullfile(deltaRoot, ['PSI ' key]);
    if ~exist(psiFolder,'dir'), mkdir(psiFolder); end

    x_phi   = PHI(mask);
    y_delta = DELTA(mask);
    z_ar    = AR(mask);

    % simple domain guard for plotting
    good = isfinite(x_phi) & isfinite(y_delta) & isfinite(z_ar);
    x_phi   = x_phi(good);
    y_delta = y_delta(good);
    z_ar    = z_ar(good);
    if isempty(x_phi), continue; end

    % --- plot 3D scatter (Z=area_ratio, X=phi, Y=ΔJcap) ---
    h = figure('Visible','on'); hold on; grid on; view(3);
    scatter3(x_phi, y_delta, z_ar, 28, 'filled');
    xlabel('\phi', 'FontSize',40, 'FontWeight','bold');
    ylabel('\Delta J_{prot} = |J_2 - J_1|', 'FontSize',40, 'FontWeight','bold');
    zlabel('area ratio @ F_{min}', 'FontSize',40, 'FontWeight','bold');
    title(sprintf('Area ratio: Z=area\\_ratio, X=\\phi, Y=\\Delta J_{prot} (\\psi = %s)', key), ...
          'FontSize',14, 'FontWeight','bold');

    savefig(h, fullfile(psiFolder, sprintf('AreaRatio_vs_phi_vs_deltaJcap_PSI_%s.fig', key)));
    saveas(h, fullfile(psiFolder, sprintf('AreaRatio_vs_phi_vs_deltaJcap_PSI_%s.png', key)));
    close(h);

    % --- save points CSV for this ψ ---
    T = table(x_phi, y_delta, z_ar, 'VariableNames', {'phi','deltaJcap','area_ratio'});
    writetable(T, fullfile(psiFolder, sprintf('points_PSI_%s.csv', key)));
end

fprintf('ΔJcap area-ratio figures written to:\n  %s\n', deltaRoot);
end


%% =======================================================================
function build_symmetricJcap_scan_AR(structPath, output_path)
% BUILD_SYMMETRICJCAP_SCAN_AR
% One global figure of all symmetric Jcap sets (Ja == Jb):
%   Z = percent deviation from analytic area
%   X = phi
%   Y = Jcap value
% Styling is the same as postproc_area_ratio_selected_pair_struct.

    if ~exist(structPath,'file')
        error('Struct file not found: %s', structPath);
    end

    % ---- load struct ----
    S = load(structPath);
    fn = fieldnames(S);
    if numel(fn) ~= 1 || ~isstruct(S.(fn{1}))
        error('File %s must contain ONE struct variable.', structPath);
    end
    S = S.(fn{1});

    hasAR_Fmin = isfield(S, 'area_ratio_at_Fmin');
    hasAR      = isfield(S, 'area_ratio');
    if ~hasAR_Fmin && ~hasAR
        error('Struct has neither "area_ratio_at_Fmin" nor "area_ratio".');
    end

    % ---- collect ALL symmetric data (Ja == Jb) ----
    phiSym   = [];
    JcapSym  = [];
    Zsym_pct = [];   % Z in percent

    for k = 1:numel(S)
        jcaps = unique(S(k).Jcap(:));
        jcaps = sort(jcaps);

        % only symmetric sets: one unique Jcap value
        if numel(jcaps) ~= 1
            continue;
        end
        Jval = jcaps(1);

        phi_k = S(k).phi(:);
        if hasAR_Fmin
            Zk = S(k).area_ratio_at_Fmin(:);
        else
            Zk = S(k).area_ratio(:);
        end

        good = isfinite(phi_k) & isfinite(Zk);
        if ~any(good), continue; end

        phi_k   = phi_k(good);
        Zk      = Zk(good);
        Zk_pct  = 100 .* Zk;   % <-- percent

        phiSym   = [phiSym;   phi_k];
        JcapSym  = [JcapSym;  repmat(Jval, numel(phi_k), 1)];
        Zsym_pct = [Zsym_pct; Zk_pct];
    end

    if isempty(phiSym)
        warning('No symmetric Jcap data (Ja == Jb) found in %s', structPath);
        return;
    end

    % ---- output folder ----
    symFolder = fullfile(output_path, 'symmetricJcaps_percent');
    if ~exist(symFolder,'dir'), mkdir(symFolder); end

    % ---- figure: same style as pair plots ----
    hFig = figure('Visible','on'); hold on; grid on; view(3);

    % ===== Half-transparent interpolated surface over (phi, Jcap) ======
    phiFit = phiSym(:);
    Jfit   = JcapSym(:);
    Zfit   = Zsym_pct(:);

    maskValid = isfinite(phiFit) & isfinite(Jfit) & isfinite(Zfit);
    phiFit = phiFit(maskValid);
    Jfit   = Jfit(maskValid);
    Zfit   = Zfit(maskValid);

    % collapse duplicates in (phi,Jcap)
    [XYu, uniqIdx] = unique([phiFit Jfit], 'rows');
    phiFit = XYu(:,1);
    Jfit   = XYu(:,2);
    Zfit   = Zfit(uniqIdx);

    if numel(Zfit) >= 3
        nGrid = 120;
        phiL = min(phiFit); phiU = max(phiFit);
        JL   = min(Jfit);   JU   = max(Jfit);

        [PhiG, JG] = meshgrid(linspace(phiL,phiU,nGrid), ...
                              linspace(JL,JU,  nGrid));

        method = 'v4';  % biharmonic spline interpolation
        ZG = griddata(phiFit, Jfit, Zfit, PhiG, JG, method);

        % keep only inside convex hull
        K  = convhull(phiFit, Jfit);
        in = inpolygon(PhiG, JG, phiFit(K), Jfit(K));
        ZG(~in) = NaN;

        set(gcf,'Renderer','opengl');
        sH = surf(PhiG, JG, ZG, ...
            'EdgeColor','none', ...
            'FaceAlpha',0.50, ...
            'DisplayName','Interpolated surface');
        shading interp;
        material dull; lighting gouraud; camlight headlight;
    end

    % ===== Color scale bar in percent, same logic as pair plots =====
if exist('sH','var') && isgraphics(sH)
    zSurfMin = min(ZG(:), [], 'omitnan');
    zSurfMax = max(ZG(:), [], 'omitnan');
    if isfinite(zSurfMin) && isfinite(zSurfMax) && zSurfMin < zSurfMax
        caxis([zSurfMin zSurfMax]);
    end

    cb = colorbar('Location','eastoutside');

    % ticks on LEFT side of the bar
    cb.YAxisLocation = 'left';

    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 22;

    % ONLY bottom + top ticks
    cb.Ticks = [zSurfMin, zSurfMax];

    % 2 decimals, percent formatting
    dec = 2;
    fmt = sprintf('%%.%df%%%%', dec);
    cb.TickLabels = {sprintf(fmt, zSurfMin), sprintf(fmt, zSurfMax)};

    set(gcf,'Renderer','opengl');
end


    % ===== Points (percent Z) – black diamonds, like AB =====
    scatter3(phiSym, JcapSym, Zsym_pct, ...
             42, 'k', 'd', 'filled', ...
             'DisplayName','symmetric J_{prot} sets');

annotation(gcf,'textbox',[0.02 0.35 0.10 0.30], ...   % [x y w h] (normalized)
    'String', { ...
        '$\frac{\left(\mathrm{micro\ area}-\mathrm{macro\ area}\right)}{\mathrm{macro\ area}}$', ...
        '', ...
        '$[\%]$' ...
    }, ...
    'Interpreter','latex', ...
    'FontSize',40, ...
    'FontWeight','bold', ...
    'Color','red', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');


    % ===== Axes, labels, title – SAME STYLE =====
    xlabel('\phi', 'FontSize',40, 'Color','red', 'FontWeight','bold');
    ylabel('$J_{\mathrm{prot}}\;\left[\frac{1}{nm}\right]$', ...
           'Interpreter','latex', ...
           'FontSize',40, ...
           'Color','red', ...
           'FontWeight','bold');

    zlabel('Area deviation (%)', ...
       'Interpreter','tex', ...
       'FontSize',40, 'Color','red', 'FontWeight','bold');

title(sprintf(['Area deviation ($\\frac{\\left(\\mathrm{micro\\ area}-\\mathrm{macro\\ area}\\right)}{\\mathrm{macro\\ area}}\\;[\\%%]$) ' ...
               ' - symmetric $J_{\\mathrm{prot}}$ pairs']), ...
      'Interpreter','latex', ...
      'FontSize',30, 'Color','red', 'FontWeight','bold');


    legend('Location','best');

    % ===== Save figure + data =====
    savefig(hFig, fullfile(symFolder, 'AreaRatio_phi_Jcap_symmetric.fig'));
    saveas(hFig, fullfile(symFolder, 'AreaRatio_phi_Jcap_symmetric.png'));
    close(hFig);

    symData = struct( ...
        'phi',            phiSym(:), ...
        'Jcap',           JcapSym(:), ...
        'area_ratio_pct', Zsym_pct(:) );

    save(fullfile(symFolder, 'symmetricJcap_points.mat'), 'symData','-v7');

    T = table(symData.phi, symData.Jcap, symData.area_ratio_pct, ...
        'VariableNames', {'phi','Jcap','area_ratio_pct'});
    writetable(T, fullfile(symFolder, 'symmetricJcap_points.csv'));

    fprintf('Symmetric Jcap percent-deviation plot written to:\n  %s\n', symFolder);
end

