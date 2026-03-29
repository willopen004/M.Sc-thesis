close all;  
clear;      
clc;
set(groot, ...
    'defaultAxesFontSize', 28, ...
    'defaultLegendFontSize', 12);   % keep legend smaller
output_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon_cylinder\V2\my_post_proc_figures\bistabilty\';
original_data_path = 'D:\yiftach_OneDrive\OneDrive - Technion\kaiser\my_output_directory_parallel\half_hexagon_cylinder\V2/families_table.mat';

main_run_function (output_path, original_data_path);

function main_run_function (output_path, original_data_path)
% -------------------------------------------------------------------------
% POST_PROCESSING_Bi_stability.M   (rev-1)
% -------------------------------------------------------------------------
% Description
% -----------
% Readable refactor of the original post-processing script with:
%   • Descriptive variable names
%   • Clear section headers and inline comments
%   • **No** change to numerical results or control flow
% -------------------------------------------------------------------------
%% ───────────────────── User-Toggles & Constants ────────────────────────
REPROCESS_TABLE   = true;          % regenerate slim table from ALL_DATA.mat

THERMAL_ENERGY    = 1;             % kBT (∆ε in the original script)
PROT_BASE_AREA    = pi * 5^2;      % nm²  – projected area of ONE protein (a₀)

ENERGY_SCALING    = 1e-10;         % converts SE output to kT
BENDING_RIGIDITY  = 20;            % κ  [kT]


%% ─────────────────────── Prepare/Load Data ─────────────────────────────
if REPROCESS_TABLE
    load(original_data_path);
    
    familiesTable.PSI_honey = 1 - familiesTable.PSI_honey;%#ok<LOAD>

    % familiesTable = familiesTable(familiesTable.PHI >= 0.05, :);

    % % Sanity check: P2 and P3 must share the same curvature
    % if any(familiesTable.Curv_P1 - familiesTable.Curv_P2)
    %     error('P2 and P3 are expected to have identical curvature');
    % end

    % % ------ Build a slim table with only the required information -------
    slim.Jcap        = [familiesTable.Curv_P1  familiesTable.Curv_P2];            % [Ja , Jb]
    % protBaseArea     = [familiesTable.protein_base_area_A  familiesTable.protein_base_area_C];        % areas of those caps
    % 
    % [~, idxSmallCap] = min(slim.Jcap, [], 2);                                 % 1 → Ja is smaller, 2 → Jb is smaller
    % 
    % slim.totalProtArea = sum(protBaseArea, 2);                                % A_p = A_a + A_b
    % slim.totalArea     = resultTable.membrane_area + slim.totalProtArea;         % A_tot = A_mem + A_p

    slim.phi = familiesTable.PHI;                          % packing fraction  φ = A_p / A_tot
    slim.f   = BENDING_RIGIDITY * ENERGY_SCALING * familiesTable.plane_bending_energy; % bending-energy density  f            
    slim.J   = familiesTable.Avg_curv;                                           % average curvature  J

    % ψ : fraction of the SMALLER cap in the mixed pair
    % slim.psi = protBaseArea(sub2ind(size(protBaseArea), (1:size(protBaseArea,1))', idxSmallCap));
    slim.psi = familiesTable.PSI_honey;

    % For symmetric caps ψ is undefined → mark as NaN
    isSymmetricPair = diff(slim.Jcap,1,2) == 0;
    slim.psi(isSymmetricPair) = NaN;

    save slim_data.mat slim                                               %#ok<SAVE>
else
    load slim_data.mat slim                                              %#ok<LOAD>
end

%% ─────────────────── Identify Unique Curvature Pairs ───────────────────
% Sort Ja/Jb per row so [Ja , Jb] is always ascending; preserve entry order
[~, firstIdx, pairIdx] = unique(sort(slim.Jcap')', 'rows', 'stable');
uniquePairs    = slim.Jcap(firstIdx, :);  % all distinct [Ja , Jb]

idxSymPairs    = find(abs(diff(uniquePairs,1,2)) == 0); % Ja == Jb
idxMixedPairs  = find(abs(diff(uniquePairs,1,2)) ~= 0); % Ja ≠ Jb

mixedPairs     = uniquePairs(idxMixedPairs, :);

%% ─────────────── Iterate over ALL Mixed Pairs ──────────────────────────
for idxPairToPlot = 1:size(mixedPairs,1)

    fprintf('\nWorking on caps: %g   %g\n', mixedPairs(idxPairToPlot, :));

    Ja = min(mixedPairs(idxPairToPlot, :));
    Jb = max(mixedPairs(idxPairToPlot, :));

    % Row selection ------------------------------------------------------
    isMixed = pairIdx == idxMixedPairs(idxPairToPlot);
    isJaJa  = pairIdx == idxSymPairs(uniquePairs(idxSymPairs,1) == Ja);
    isJbJb  = pairIdx == idxSymPairs(uniquePairs(idxSymPairs,1) == Jb);
    if ~any(isJaJa) || ~any(isJbJb)
        error('Missing symmetric data for this pair');
    end

    % Assemble vectors ---------------------------------------------------
    J   = [slim.J(isMixed);  slim.J(isJaJa);   slim.J(isJbJb)];
    F   = [slim.f(isMixed);  slim.f(isJaJa);   slim.f(isJbJb)];
    phi = [slim.phi(isMixed); slim.phi(isJaJa); slim.phi(isJbJb)];
    psi = [slim.psi(isMixed); ones(sum(isJaJa),1); zeros(sum(isJbJb),1)];

    % Common sampling parameters ----------------------------------------
    NUM_J_SAMPLES   = 200;
    NUM_PSI_SAMPLES = 200;
    NUM_PHI_SLICES  = 10;

    phiGrid = linspace(min(phi), max(phi), NUM_PHI_SLICES).';
    Jgrid   = linspace(min(J),   max(J),   NUM_J_SAMPLES);   % <-- ONE grid!

    % Pre-allocate storage for subplot-3 data
    FoptSlices  = zeros(NUM_PHI_SLICES, NUM_J_SAMPLES);
    psiOptSlices = zeros(NUM_PHI_SLICES, NUM_J_SAMPLES);     % if you need ψopt

    % ===================================================================
    % φ-LOOP
    % ===================================================================
    for p = NUM_PHI_SLICES:-1:1
        phiSlice = phiGrid(p);

        inSlice = phi >= phiSlice*(1-0.5/NUM_PHI_SLICES) & ...
                  phi <= phiSlice*(1+0.5/NUM_PHI_SLICES);

        fprintf('φ-slice %.3g → keeping %d / %d points\n', ...
        phiSlice, sum(inSlice), numel(inSlice));

        Jdata   = J(inSlice);
        psidata = psi(inSlice);
        Fdata   = F(inSlice);
        % ----- Filter current φ-slice ------------------------------------------
        valid    = ~isnan(Jdata) & ~isnan(psidata) & ~isnan(Fdata);
        Jdata    = Jdata(valid);
        psidata  = psidata(valid);
        Fdata    = Fdata(valid);
        
        if numel(Jdata) < 6                    % need ≥6 finite points for poly22
            fprintf('φ-slice %.4g   skipped (%d valid points)\n', ...
                    phiSlice, numel(Jdata));
            continue                           % jump to next φ-slice
        end
        
        % ---------- now it is safe to build the meshes --------------------------

        % Surface fit f(J,ψ)
        Jmesh   = linspace(min(Jdata), max(Jdata), NUM_J_SAMPLES);
        psimesh = linspace(0,1,NUM_PSI_SAMPLES);
        [JJ,PP] = meshgrid(Jmesh, psimesh);
        

        surfFit = fit([Jdata, psidata], Fdata, 'poly22');
        deltaF  = min(Fdata) - min(surfFit(JJ,PP),[],'all');

        bendE = @(Jv,psv) surfFit(Jv,psv) + deltaF;
        entE  = @(psv) (phiSlice/PROT_BASE_AREA) .* ...
                       (psv.*log(psv) + (1-psv).*log(1-psv));
        totF  = @(Jv,psv) bendE(Jv,psv) + entE(psv);

        % ---------- minimise along ψ on the COMMON Jgrid ---------------
        for k = 1:NUM_J_SAMPLES
            [psiOptSlices(p,k), FoptSlices(p,k)] = ...
                fminbnd(@(y) totF(Jgrid(k),y), 0, 1);
        end

        % ---------- individual slice figure (unchanged) ----------------
        figure('Name',sprintf('φ = %.3f   Ja=%.3g  Jb=%.3g',phiSlice,Ja,Jb));

        subplot(3,1,1)
        scatter3(Jdata,psidata,Fdata,'filled'); hold on
        surf(JJ,PP,bendE(JJ,PP),'EdgeColor','none','FaceColor','r', ...
             'FaceAlpha',0.5); hold off
        view(30,25); xlabel('J'); ylabel('\psi'); zlabel('f');
        title(sprintf('\\phi = %.3f',phiSlice))

        subplot(3,1,2)
        surf(JJ,PP,totF(JJ,PP),'EdgeColor','none'); view(30,25);
        xlabel('J'); ylabel('\psi'); zlabel('f_{total}');

        subplot(3,1,3)
        plot(Jgrid,FoptSlices(p,:),'LineWidth',4.5);
        xlabel('J [nm^{-1}]'); ylabel('f(\\psi_{opt}) [kT nm^{-2}]'); grid on;

        title(sprintf('Ja = %.3g   Jb = %.3g   \\phi = %.3f', Ja, Jb, phiSlice));

        % pause;
            % (optional) save the figure
        pairFolder = fullfile(output_path, ...
                              sprintf('Jcap_combination[%g,%g]',Ja,Jb));
        if ~exist(pairFolder,'dir'); mkdir(pairFolder); end
        % savefig(gcf, fullfile(pairFolder,sprintf('Phi_Slice - %f.fig', phiSlice)));
        close(gcf)   % comment out if you want to keep each slice figure
    end

    phiLabels = arrayfun(@(v) sprintf('\\phi = %.3g', v), phiGrid, 'uni', 0);

% ===================================================================
% Combined 3-D plot  Z = fopt , X = J , Y = φ
% ===================================================================
figure('Name',sprintf('f(ψ_{opt}) vs J & φ   Ja=%.3g  Jb=%.3g',Ja,Jb));
hold on
L = gobjects(NUM_PHI_SLICES,1);
for p = 1:NUM_PHI_SLICES
    L(p) = plot3(Jgrid, phiGrid(p)*ones(1,NUM_J_SAMPLES), FoptSlices(p,:), ...
                 'LineWidth',4.5, ...
                 'DisplayName',phiLabels{p});
end
hold off
xlabel('$J\;\left[\frac{1}{nm}\right]$', ...
           'Interpreter','latex', ...
           'FontSize',40, ...
           'Color','red', ...
           'FontWeight','bold');
ylabel('\phi', 'FontSize',40, 'Color','red');
zlabel('$f(\psi_{\mathrm{opt}})\;\left[\frac{k_{\mathrm{B}}T}{\mathrm{nm}^{2}}\right]$', ...
       'Interpreter','latex', ...
       'FontSize',40, ...
       'Color','red', ...
       'FontWeight','bold');

title(sprintf('f(ψ_{opt}) vs J vs \\phi [%.3g, %.3g]',Ja,Jb),'Color','red');
grid on; view(45,30); box on;
legend(L,phiLabels,'Location','best','Interpreter','tex');

pairFolder = fullfile(output_path,sprintf('Jcap_combination[%g,%g]',Ja,Jb));
if ~exist(pairFolder,'dir'); mkdir(pairFolder); end
savefig(gcf, fullfile(pairFolder,'Fopt_vs_J_and_phi.fig'));
close(gcf);


% ──────────────────────────────────────────────────────────────────────
% New 3-D plot:  Z = ψopt    |    X = J    |    Y = φ
% ──────────────────────────────────────────────────────────────────────
figure('Name',sprintf('ψ_{opt} vs J & φ   Ja=%.3g  Jb=%.3g',Ja,Jb));
hold on
L2 = gobjects(NUM_PHI_SLICES,1);
for p = 1:NUM_PHI_SLICES
    L2(p) = plot3(Jgrid, phiGrid(p)*ones(1,NUM_J_SAMPLES), psiOptSlices(p,:), ...
                  'LineWidth',4.5, ...
                  'DisplayName',phiLabels{p});
end
hold off
xlabel('$J\;\left[\frac{1}{nm}\right]$', ...
           'Interpreter','latex', ...
           'FontSize',40, ...
           'Color','red', ...
           'FontWeight','bold');
ylabel('\phi', 'FontSize',40, 'Color','red');
zlabel('ψ_{opt}', 'FontSize',40, 'Color','red');
title(sprintf('ψ optimal vs J vs \\phi [%.3g, %.3g]', Ja, Jb), ...
        'Color','red', ...
        'FontSize',30, ...
        'FontWeight','bold');

grid on; view(45,30); box on;
legend(L2,phiLabels,'Location','best','Interpreter','tex');

pairFolder = fullfile(output_path,sprintf('Jcap_combination[%g,%g]',Ja,Jb));
if ~exist(pairFolder,'dir'); mkdir(pairFolder); end
savefig(gcf, fullfile(pairFolder,'PSIopt_vs_J_and_phi.fig'));
close(gcf);


end

end

