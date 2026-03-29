close all;
clear all;
clc;

%% this script is for running with larger PHI of all proteins

template_dir = 'E:\Yiftach\OneDrive - Technion\work_folder\clean scripts\hexagon\data collection\template for surface evolver';  % Specify the directory where your template file is located
output_dir = 'E:\Yiftach\OneDrive - Technion\kaiser\my_output_directory_parallel\half_triangle\V1\';  % Specify the directory where you want to save the results

Min_Dist_var = 25;
Max_Dist_var = 45;
Prot_rad_var = 5;
AVG_CURVATURE_var = 0.001;
make_plots = false; % display the Dist distribution
num_PHI_values = 15; % number of Dist values running
sim_only = true;
RATIO_var = [1, 11, 13 , 15];
% NUMBER_OF_PROTEINS_var = [0, 10, 42, 106, 204];
% original Jcaps:
% Jcap1 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15]; 
% Jcap2 = [0.001, 0.005, 0.016, 0.031, 0.046, 0.061, 0.076, 0.09, 0.105, 0.12, 0.135, 0.15];
% Jcap1 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];
% Jcap2 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];
Jcap1 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];
Jcap2 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];
curvature_value_range1 = NaN;

% result = run_single_evolver_iteration(204, 0.00001, 30, 0.15, 0.001, template_dir, output_dir);

% run_model (Dist_var, P_var, AVG_CURVATURE_var, NUMBER_OF_PROTEINS_var);

stats = run_model(Min_Dist_var, Max_Dist_var, Prot_rad_var, AVG_CURVATURE_var, RATIO_var, num_PHI_values, make_plots);

check_curvature_combinations( ...
    Min_Dist_var, Max_Dist_var, Prot_rad_var, AVG_CURVATURE_var, ...
    num_PHI_values, make_plots, ...
    RATIO_var, Jcap1, Jcap2, output_dir, template_dir, Prot_rad_var, sim_only, curvature_value_range1);

function check_curvature_combinations( ...
    Min_Dist_var, Max_Dist_var, Prot_rad_var, AVG_CURVATURE_var, ...
    num_PHI_values, make_plots, ...
    RATIO_var, Jcap1, Jcap2, results_folder, template_dir, rad_P, sim_only, curvature_value_range1)

% Rules:
% • If Jcap1 == Jcap2  → ONLY ratio=1 and ONLY protein=0.
% • Else               → ratios in RATIO_var using mapping:
%       11 → [42]
%       13 → [106, 204]
%       15 → [106, 10, 204]
% Stats:
% • Recompute stats **per ratio** via run_model(..., ratio_value, ...).
% Skips if any of these exist:
%   1) results_folder/RATIO_<ratio>/model_<P>_proteins/results_<P>_...
%   2) results_folder/RATIO_<ratio>/model_<P>_proteins/ratio_<r>__results_... (future)
%   3) results_folder/model_<P>_proteins/results_<P>_... (legacy)

    if nargin < 14 || isempty(sim_only)
        sim_only = false;   % default: process all pairs (as before)
    end
    validateattributes(sim_only, {'logical','numeric'}, {'scalar'}, mfilename, 'sim_only');
    sim_only = logical(sim_only);

    if isempty(Jcap1) || isempty(Jcap2)
        error('Jcap1 and Jcap2 must not be empty');
    end
    if isempty(RATIO_var)
        warning('RATIO_var is empty. Nothing to run.');
        return;
    end

    tol = 1e-12;

    for j = 1:numel(Jcap1)
        for k = 1:numel(Jcap2)
            J1 = Jcap1(j);  J2 = Jcap2(k);
            J1_str = sprintf('%.3f', J1);
            J2_str = sprintf('%.3f', J2);

            % -------- optional filter for "similar-only" mode --------
            if sim_only && abs(J1 - J2) > tol
                fprintf('• sim_only=TRUE → skip unequal pair [%s,%s].\n', J1_str, J2_str);
                continue;
            end

            % -------- pair-level routing --------
            if abs(J1 - J2) <= tol
                ratios_here    = 1;   % ONLY ratio=1
                fixed_proteins = 0;   % ONLY protein=0
                fprintf('• Pair [%s,%s] equal → restrict to ratio=1, protein=0.\n', J1_str, J2_str);
            else
                ratios_here    = setdiff(RATIO_var, 1, 'stable'); % drop 1 for unequal pairs
                fixed_proteins = [];
            end

            for r = 1:numel(ratios_here)
                ratio_value = ratios_here(r);
                ratio_str   = toChar(ratio_value);

                % Proteins for this ratio
                if ~isempty(fixed_proteins)
                    proteins_to_run = fixed_proteins;   % equal-pair case
                else
                    proteins_to_run = proteins_for_ratio(ratio_value);
                end
                if isempty(proteins_to_run)
                    fprintf('• Ratio=%s → no proteins mapped; skipping pair [%s,%s].\n', ...
                            ratio_str, J1_str, J2_str);
                    continue;
                end

                % ---- Check existing target folders FIRST (per protein) ----
                ratio_root_new = fullfile(results_folder, sprintf('RATIO_%s', ratio_str));

                need_runs_mask = false(size(proteins_to_run));
                for idx = 1:numel(proteins_to_run)
                    P = proteins_to_run(idx);

                    % (1) Current scheme
                    model_dir_new        = fullfile(ratio_root_new, sprintf('model_%d_proteins', P));
                    sub_unprefixed       = sprintf('results_%d_proteins_curvatures_[%s, %s]', P, J1_str, J2_str);
                    full_new_unprefixed  = fullfile(model_dir_new, sub_unprefixed);

                    % (2) Future scheme with ratio prefix
                    sub_prefixed         = sprintf('ratio_%s__results_%d_proteins_curvatures_[%s, %s]', ...
                                                   ratio_str, P, J1_str, J2_str);
                    full_new_prefixed    = fullfile(model_dir_new, sub_prefixed);

                    % (3) Legacy scheme (no RATIO layer)
                    model_dir_legacy     = fullfile(results_folder, sprintf('model_%d_proteins', P));
                    full_legacy          = fullfile(model_dir_legacy, sub_unprefixed);

                    if isfolder(full_new_unprefixed) || isfolder(full_new_prefixed) || isfolder(full_legacy)
                        fprintf('    ✔ SKIP: found existing folder for RATIO:%s P=%d [%s,%s].\n', ...
                                ratio_str, P, J1_str, J2_str);
                        need_runs_mask(idx) = false;
                    else
                        fprintf('    ↻ WILL RUN: RATIO:%s P=%d [%s,%s].\n', ...
                                ratio_str, P, J1_str, J2_str);
                        need_runs_mask(idx) = true;
                    end
                end

                if ~any(need_runs_mask)
                    % Nothing to do for this (pair, ratio)
                    continue;
                end

                % ---- RECOMPUTE STATS *PER RATIO* (exactly once here) ----
                % run_model takes the ratio as the selector that drives PSI etc.
                base_stats = run_model( ...
                    Min_Dist_var, Max_Dist_var, Prot_rad_var, AVG_CURVATURE_var, ...
                    ratio_value, num_PHI_values, make_plots);

                % Build a stats array for this ratio, keyed by the proteins we’ll run
                stats_curr = repmat(base_stats(1), 1, numel(proteins_to_run));
                for ii = 1:numel(proteins_to_run)
                    stats_curr(ii).PROTEINS = proteins_to_run(ii);  % map key for iterative_evolver_run
                    stats_curr(ii).RATIO    = ratio_value;          % keep the ratio with the stats
                end

                % ---- Launch only the ones that need runs ----
                todo = proteins_to_run(need_runs_mask);
                fprintf('• Ratio=%s → RUN proteins [%s] | pair [%s, %s]\n', ...
                        ratio_str, strjoin(string(todo), ','), J1_str, J2_str);

                % Ensure parent model directories exist
                for P = todo
                    mdl = fullfile(ratio_root_new, sprintf('model_%d_proteins', P));
                    if ~isfolder(mdl), try, mkdir(mdl); catch, end, end
                end

                % Pass ratio-root so outputs land under RATIO_<ratio>/
                for P = todo
                    activate_and_monitor_parallel( ...
                        template_dir, ratio_root_new, J1, J2, P, ...
                        rad_P, stats_curr, curvature_value_range1, ratio_value);
                end

                % Explicitly clear per-ratio stats to avoid accidental reuse
                clear stats_curr base_stats

            end % ratio loop
        end % Jcap2
    end % Jcap1
end

% ---------- helpers ----------
function vec = proteins_for_ratio(ratio_value)
    switch ratio_value
        case 11
            vec = [42];
        case 13
            vec = [106, 204];
        case 15
            vec = [106, 10, 204];
        otherwise
            vec = []; % ratio 1 is handled only in equal-pair rule
    end
end

function c = toChar(x)
    if isstring(x) || ischar(x)
        c = char(x);
    elseif isnumeric(x)
        c = num2str(x);
    else
        try, c = char(string(x));
        catch, c = '<unk>';
        end
    end
end




function activate_and_monitor_parallel(template_dir, output_dir, P1_curv_values_var, P2_P3_curv_values_var1, num_proteins, rad_P, stats, curvature_value_range, RATIO_VAR)
    % Ensure the Parallel Computing Toolbox is available
    if isempty(ver('parallel'))
        error('Parallel Computing Toolbox is not installed.');
    end

    % Check if a parallel pool is already running
    pool = gcp('nocreate');
    if isempty(pool)
        % Start a parallel pool with max available workers
        pool = parpool('local');
    else
        fprintf('Parallel pool already running with %d workers.\n', pool.NumWorkers);
    end

    % Log file for monitoring
    log_file = fullfile(pwd, 'parallel_execution_log.txt');

    % Ensure the directory for log file exists
    log_folder = fileparts(log_file);
    if ~exist(log_folder, 'dir')
        mkdir(log_folder);
    end

    if exist(log_file, 'file')
        delete(log_file); % Clear previous log
    end

    % Start monitoring system performance
    fprintf('Starting parallel execution monitoring...\n');
    diary(log_file); % Redirect output to log file

    % Time execution
    tic;
    try
        % Call the main function
        iterative_evolver_run(template_dir, output_dir, P1_curv_values_var, P2_P3_curv_values_var1, num_proteins, rad_P, stats, curvature_value_range, RATIO_VAR)
    catch ME
        fprintf('Error in execution: %s\n', ME.message);
        
        % Display the function and line number where the error occurred
        for k = 1:length(ME.stack)
            fprintf('Error in function: %s (Line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
    end
    elapsed_time = toc;

    % Stop monitoring
    diary off;

    % Display final results
    fprintf('Parallel execution completed in %.2f seconds.\n', elapsed_time);
    fprintf('Check log file: %s\n', log_file);
end


function iterative_evolver_run(template_dir, output_dir, ...
    P1_curv_values, P2_P3_curv_values, ...
    num_proteins_vec, rad_P, stats, curvature_value_range, RATIO_VAR)
% ITERATIVE_EVOLVER_RUN  – launch Evolver jobs and save .MAT results only.
% All curvature samples *and* the distance grid come straight from `stats`.

% -----------------------------------------------------------------------
% ------------------------------------------------------------------
% Build map <protein‑count → index in stats>  (skip empties / dups)
% ------------------------------------------------------------------
valid = arrayfun(@(s) isnumeric(s.PROTEINS) && isscalar(s.PROTEINS) ...
                              && ~isempty(s.PROTEINS), stats);

allKeys = arrayfun(@(s) s.PROTEINS, stats(valid));   % numeric row vector
allVals = find(valid);                               % matching indices

% remove duplicates, keep first occurrence ("stable" order)
[uniqKeys, ia] = unique(allKeys, 'stable');
uniqVals       = allVals(ia);

statsMap = containers.Map(num2cell(uniqKeys), num2cell(uniqVals));
if ~exist(output_dir,'dir'), mkdir(output_dir); end
% -----------------------------------------------------------------------

for N = num_proteins_vec
    if ~statsMap.isKey(N), warning('No stats for N=%d – skip.',N); continue; end
    sIdx     = statsMap(N);
    distGrid = stats(sIdx).Dist_grid;      % distance vector for this N
    nDist    = numel(distGrid);

    protDir = fullfile(output_dir,sprintf('model_%d_proteins',N));
    if ~exist(protDir,'dir'), mkdir(protDir); end

    for J1 = P1_curv_values
        for J23 = P2_P3_curv_values

            fprintf('Prot %d  P1=%.3f  P23=%.3f  (#dists=%d)\n',N,J1,J23,nDist);
            pairDir = fullfile(protDir, ...
                sprintf('results_%d_proteins_curvatures_[%.3f, %.3f]', N, J1, J23));
            if ~exist(pairDir,'dir'), mkdir(pairDir); end

            % ------------- build curvature samples -------------------
            resCurv = analyzeProteinCurvature(stats,N,J1,J23); % no distRanges
            CurvMat = resCurv.J_range_values;                  % Ndists × 20

            % ------------- distance loop -----------------------------
            for dIdx = 1:nDist
                Dist_val = distGrid(dIdx);
                
                % ▼── choose where CurvRow comes from
                if ~isnan(curvature_value_range(1))
                    CurvRow = curvature_value_range;          % external override
                else
                    CurvRow = unique(CurvMat(dIdx,:),'stable'); % default row from analyze*
                end

                nCurv    = numel(CurvRow);

                fprintf('[ROW] Prot %d | Dist %.2f | PHI %.6f | Curv %.10f → %.10f\n',...
                        N,Dist_val,resCurv.PHI_values(dIdx),CurvRow(1),CurvRow(end));

                                % ---- launch parfeval jobs ----------------------------
                pool = gcp();
                futs = parallel.FevalFuture.empty(nCurv,0);
                startUTC = NaT(1,nCurv,"TimeZone","UTC");

                PHI_now  = resCurv.PHI_values(dIdx);
                Jint_cent= resCurv.Jint_values(dIdx);   % centre for this φ row
                PSI_val  = resCurv.PSI;                % constant for this N
                Jint_theory = PHI_now * (PSI_val * J1 + (1-PSI_val) * J23);

                for k = 1:nCurv
                    CurvVal = CurvRow(k);

                    % ---------- DEBUG line ---------------------------
                    fprintf(['[DEBUG] Prot %d | Dist %.2f | Curv %.6f | PHI %.6f | PSI %.4f | ' ...
                             'Jint %.6f | P1 %.3f | P23 %.3f  =>  Jint theory %.8f = %.6f * (%.4f * %.3f + (1-%.4f) * %.3f)\n'], ...
                            N,           Dist_val, CurvVal, PHI_now, PSI_val, ...
                            Jint_cent,   J1,       J23,     Jint_theory, ...
                            PHI_now,     PSI_val,  J1,      PSI_val,     J23);

                    % ---------- launch worker ------------------------
                    futs(k) = parfeval(pool,@run_single_evolver_iteration,1, ...
                               N, CurvVal, Dist_val, J1, J23, RATIO_VAR, PSI_val, template_dir, pairDir);
                end


                % -- timeout watchdog --------------------------------
                timeoutSec = 300; poll = 0.25;
                resStruct  = repmat(buildNaNStruct(rad_P,J1,J23,NaN,Dist_val,N),1,nCurv);
                
                % --- inject the curvature value for each row --------------------------
                for kk = 1:nCurv
                    resStruct(kk).Avg_curv = CurvRow(kk);   % <‑‑ keep even on timeout
                end

                hasStarted = false(1,nCurv);
                isClosed   = false(1,nCurv);
                
                while ~all(isClosed)
                    % A) try to fetch a finished future (times out after poll seconds)
                    idxDone = [];
                    resultData = [];
                    try
                        [idxDone,resultData] = fetchNext(futs,poll);
                    catch
                        idxDone = find(strcmp({futs.State},"finished") & ...
                                       ~cellfun(@isempty,{futs.Error}) & ~isClosed,1);
                    end
                
                    % B) mark newly‑running jobs (start timer)
                    for k = find(~hasStarted)
                        if futs(k).State=="running"
                            hasStarted(k) = true;
                        end
                    end
                
                    % C) process finished / errored / cancelled jobs
                    if ~isempty(idxDone) && ~isClosed(idxDone)
                        runSec = seconds(futs(idxDone).RunningDuration);
                
                        if isempty(futs(idxDone).Error) && runSec < timeoutSec
                            resStruct(idxDone) = resultData;
                            fprintf('Protein %d → DistIdx %d (%.2f) ✔ Curv = %.8f (J1=%.3f, J2=%.3f) finished successfully in %.2fs\n', ...
                                    N, dIdx, Dist_val, CurvRow(idxDone), J1, J23, runSec);
                        else
                            cancel(futs(idxDone));
                            resStruct(idxDone).error_message = "timeout/error";
                            if runSec >= timeoutSec
                                fprintf(2,'Protein %d → DistIdx %d (%.2f) ⏱ Curv = %.8f (J1=%.3f, J2=%.3f) exceeded %gs (%.2fs) → NaNs\n', ...
                                        N, dIdx, Dist_val, CurvRow(idxDone), J1, J23, timeoutSec, runSec);
                            else
                                e = futs(idxDone).Error;
                                fprintf(2,'Protein %d → DistIdx %d (%.2f) ⚠ Curv = %.8f (J1=%.3f, J2=%.3f) ERROR: %s → NaNs\n', ...
                                        N, dIdx, Dist_val, CurvRow(idxDone), J1, J23, e.message);
                            end
                        end
                        isClosed(idxDone) = true;
                        drawnow
                    end
                
                    % D) cancel any running job that now exceeds timeout
                    for k = find(~isClosed & hasStarted)
                        runSec = seconds(futs(k).RunningDuration);
                        if runSec >= timeoutSec
                            cancel(futs(k));
                            resStruct(k).error_message = "timeout";
                            fprintf(2,'Protein %d → DistIdx %d (%.2f) ⏱ Curv = %.8f (J1=%.3f, J2=%.3f) exceeded %gs (%.2fs) → NaNs\n', ...
                                    N, dIdx, Dist_val, CurvRow(k), J1, J23, timeoutSec, runSec);
                            isClosed(k) = true;
                        end
                    end
                
                    % E) brief pause if some slots still pending
                    if any(~hasStarted & ~isClosed), pause(poll); end
                end
                
                % -- save MAT (table!) --------------------------------
                tbl = struct2table(resStruct(:),'AsArray',true);
                fileName = sprintf('results_%d_prot_Distance_%.2f_prot_curv_[%.3f,%.3f].mat', ...
                                   N, Dist_val, J1, J23);
                save(fullfile(pairDir, fileName), 'tbl', '-v7');

            end
        end
    end
end

disp('All iterations complete – MAT tables written.');
end

function result = run_single_evolver_iteration(current_num_proteins, Curvature, Distance, P1_curv, P2_3_curv, ratio_value, PSI_val, template_dir, output_dir)
    % Initialize result
% --- initialize result fields to mirror OUTPUTLINE names exactly ---
result.radius_P1 = NaN;
result.radius_P2 = NaN;
result.radius_P3 = NaN;
result.radius_P4 = NaN;
result.radius_P5 = NaN;

result.Curv_P1_fam = NaN;
result.Curv_P2_fam = NaN;

result.Curv_P1 = P1_curv;
result.Curv_P2 = P2_3_curv;
result.Curv_P3 = P2_3_curv;
result.Curv_P4 = P2_3_curv;
result.Curv_P5 = P2_3_curv;

result.Avg_curv = Curvature;
result.distance = Distance;
result.membrane_tension = NaN;

result.total_energy = NaN;
result.total_Bending_energy = NaN;
result.total_Bending_energy_calculated = NaN;
result.plane_bending_energy = NaN;   % NOTE: was "plain_Bending_energy"
result.cap_bending_energy = NaN;

result.total_area = NaN;
result.protein_caps_area = NaN;
result.membrane_area = NaN;
result.big_sphere_area = NaN;
result.big_sphere_calc_area = NaN;

result.num_triangles = NaN;

result.protein_base_area_A = NaN;
result.protein_base_area_B = NaN;
result.protein_base_area_C = NaN;
result.protein_base_area_D = NaN;
result.protein_base_area_E = NaN;
result.total_protein_base_area = NaN;

result.num_proteins = NaN;
result.area_per_protein_anltc = NaN;
result.area_per_protein_real  = NaN;

result.P_in_model = current_num_proteins;

result.arc_AB = NaN;
result.arc_AC = NaN;
result.arc_AD = NaN;
result.arc_AE = NaN;
result.arc_BC = NaN;
result.arc_CD = NaN;
result.arc_DE = NaN;

result.analytic_total_area = NaN;
result.PSI = NaN;
result.PSI_honey = NaN;
result.PHI = NaN;
result.RATIO = NaN;

result.error_message = '';  % keep for diagnostics

    

    % Full path to the template file in the template directory
    fname_template = fullfile(template_dir, '3_protein_triangle_config_curvature_template.fe');
    % fname_run = fullfile(output_dir, sprintf('run_curvature_configuration_%d_proteins_sphere_curv_%s_Proteins_curv_[%.3f, %.3f]_distance_%d.fe',current_num_proteins,Curvature, P1_curv, P2_3_curv, Distance));  % Save the run file in the specified directory
    fname_run = fullfile(output_dir, ...
    sprintf('run_%dP_sphCurv_%0.8f_P1_%0.3f_P23_%0.3f_dist_%.4f.fe', ...
    current_num_proteins, Curvature, P1_curv, P2_3_curv, Distance));

    % >>> ADD — begin -------------------------------------------------------------
    cmdName   = strrep(getFilename(fname_run),'fe','cmd');        % cmd twin
    fname_cmd = fullfile(output_dir, cmdName);
    cleaner = onCleanup(@() safeDeleteFiles({fname_cmd}));
    % >>> ADD — end   -------------------------------------------------------------

    fid = fopen(fname_template, 'rt');
    if fid == -1
        result.error_message = 'Could not open template file';
        return;
    end
    X = fread(fid);
    fclose(fid);

    % Replace placeholders with the current variable values
    X = char(X.');
    Y = strrep(X, 'CURV_VAR', num2str(Curvature));  % Replace Ratio placeholder
    Y = strrep(Y, 'DIST_VAR', num2str(Distance));  % Replace membrane_tension placeholder
    Y = strrep(Y, 'NUM_PROTEIN_VAR', num2str(current_num_proteins));   % Replace P_2_3_CURV_VAR with current Curvature
    Y = strrep(Y, 'Curv_P1_fam_VAR', num2str(P1_curv));   % curv P1 for the family
    Y = strrep(Y, 'Curv_P2_fam_VAR', num2str(P2_3_curv));   % curv P2 for the family
    red = '4';
    blue = '1';
        if (current_num_proteins == 0)
        
        % proteins curvature
        Y = strrep(Y, 'P_A_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_C_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_D_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P1_curv));
        
        % proteins colors
        Y = strrep(Y, 'P_A_COLOR_VAR', red);
        Y = strrep(Y, 'P_B_COLOR_VAR', red);
        Y = strrep(Y, 'P_C_COLOR_VAR', red);
        Y = strrep(Y, 'P_D_COLOR_VAR', red);
        Y = strrep(Y, 'P_E_COLOR_VAR', red);
        
        elseif (current_num_proteins == 106)
        
        % proteins curvature
        Y = strrep(Y, 'P_A_CURVA_VAR', num2str(P2_3_curv));
        Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_C_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_D_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P1_curv));
        
        % proteins colors
        Y = strrep(Y, 'P_A_COLOR_VAR', blue);
        Y = strrep(Y, 'P_B_COLOR_VAR', red);
        Y = strrep(Y, 'P_C_COLOR_VAR', red);
        Y = strrep(Y, 'P_E_COLOR_VAR', red);
        Y = strrep(Y, 'P_D_COLOR_VAR', red);
        
        elseif (current_num_proteins == 10)
        
        % proteins curvature
        Y = strrep(Y, 'P_A_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P2_3_curv));
        Y = strrep(Y, 'P_C_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_D_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P1_curv));
        
        % proteins colors
        Y = strrep(Y, 'P_A_COLOR_VAR', red);
        Y = strrep(Y, 'P_B_COLOR_VAR', blue);
        Y = strrep(Y, 'P_C_COLOR_VAR', red);
        Y = strrep(Y, 'P_D_COLOR_VAR', red);
        Y = strrep(Y, 'P_E_COLOR_VAR', red);
        
        elseif (current_num_proteins == 204)
        
        % proteins curvature
        Y = strrep(Y, 'P_A_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P2_3_curv));
        Y = strrep(Y, 'P_C_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_D_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P2_3_curv));
        
        % proteins colors
        Y = strrep(Y, 'P_A_COLOR_VAR', red);
        Y = strrep(Y, 'P_B_COLOR_VAR', blue);
        Y = strrep(Y, 'P_C_COLOR_VAR', red);
        Y = strrep(Y, 'P_D_COLOR_VAR', red);
        Y = strrep(Y, 'P_E_COLOR_VAR', blue);
        
        elseif (current_num_proteins == 42)
        
        % proteins curvature
        Y = strrep(Y, 'P_A_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P1_curv));
        Y = strrep(Y, 'P_C_CURVA_VAR', num2str(P2_3_curv));
        Y = strrep(Y, 'P_D_CURVA_VAR', num2str(P2_3_curv));
        Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P1_curv));
        
        % proteins colors
        Y = strrep(Y, 'P_A_COLOR_VAR', red);
        Y = strrep(Y, 'P_B_COLOR_VAR', red);
        Y = strrep(Y, 'P_C_COLOR_VAR', blue);
        Y = strrep(Y, 'P_D_COLOR_VAR', blue);
        Y = strrep(Y, 'P_E_COLOR_VAR', red);
        end


    % Write the updated .fe file in the specified output directory
    fid2 = fopen(fname_run, 'wt');
    if fid2 == -1
        result.error_message = 'Could not write run file';
        return;
    end
    fwrite(fid2, Y);
    fclose(fid2);

    % disp(fname_run);
    % Run the Evolver file and capture the output
    [status, out_string] = runEvolverFileParallel(fname_run, 'runme();', [], [], false);

    if status == 0
        disp('Surface Evolver ran successfully.');
        % Process the output to extract the required variables (e.g., Je, area, etc.)
        TextAsCells = regexp(out_string, '\n', 'split');
        linenum = find(~cellfun(@isempty, strfind(TextAsCells, 'OUTPUTLINE')));
            if isempty(linenum)
                % ① Populate “safe‑defaults” so the caller can still rely on result.*
                % ② Preserve basic input values so later code can tell what was run
            % ① Preserve basic inputs
            result.radius_P1 = NaN;
            result.radius_P2 = NaN;
            result.radius_P3 = NaN;
            result.radius_P4 = NaN;
            result.radius_P5 = NaN;

            result.Curv_P1_fam = Curv_P1_fam;
            result.Curv_P2_fam = Curv_P2_fam;
            
            result.Curv_P1   = P1_curv;
            result.Curv_P2   = P2_3_curv;
            result.Curv_P3   = P2_3_curv;
            result.Curv_P4   = P2_3_curv;
            result.Curv_P5   = P2_3_curv;
            
            result.Avg_curv  = Curvature;
            result.distance  = Distance;
            
            % ② Everything else = NaN/defaults, with names matching OUTPUTLINE exactly
            result.membrane_tension                 = NaN;
            result.total_energy                     = NaN;
            result.total_Bending_energy             = NaN;
            result.total_Bending_energy_calculated  = NaN;
            result.plane_bending_energy             = NaN;   % NOTE: renamed from plain_Bending_energy
            result.cap_bending_energy               = NaN;
            
            result.total_area           = NaN;
            result.protein_caps_area    = NaN;
            result.membrane_area        = NaN;
            result.big_sphere_area      = NaN;
            result.big_sphere_calc_area = NaN;
            
            result.num_triangles = NaN;
            
            result.protein_base_area_A      = NaN;
            result.protein_base_area_B      = NaN;
            result.protein_base_area_C      = NaN;
            result.protein_base_area_D      = NaN;
            result.protein_base_area_E      = NaN;
            result.total_protein_base_area  = NaN;
            
            result.num_proteins            = NaN;
            result.area_per_protein_anltc  = NaN;
            result.area_per_protein_real   = NaN;
            
            result.P_in_model = current_num_proteins;
            
            result.arc_AB = NaN;
            result.arc_AC = NaN;
            result.arc_AD = NaN;
            result.arc_AE = NaN;
            result.arc_BC = NaN;
            result.arc_CD = NaN;
            result.arc_DE = NaN;
            
            result.analytic_total_area = NaN;
            result.PSI = NaN;
            result.PSI_honey = NaN;
            result.PHI = NaN;
            result.RATIO = NaN;
            
            result.error_message = 'model converged wrong (no OUTPUTLINE)';

            else

            for n = 1:numel(linenum)
                tmp = TextAsCells{linenum(n)};
                tmp = strrep(tmp, 'OUTPUTLINE ', '');
                eval(tmp);  % Update the corresponding variables
            end
    
         % Now assign the values extracted from Surface Evolver Captured from OUTPUTLINE in Evolver

            result.radius_P1 = radius_P1;
            result.radius_P2 = radius_P2;
            result.radius_P3 = radius_P3;
            result.radius_P4 = radius_P4;
            result.radius_P5 = radius_P5;

            result.Curv_P1_fam = Curv_P1_fam;
            result.Curv_P2_fam = Curv_P2_fam;
            
            result.Curv_P1 = Curv_P1;
            result.Curv_P2 = Curv_P2;
            result.Curv_P3 = Curv_P3;
            result.Curv_P4 = Curv_P4;
            result.Curv_P5 = Curv_P5;
            
            result.Avg_curv = Avg_curv;
            result.distance = distance;
            result.membrane_tension = membrane_tension;
            
            result.total_energy = total_energy;
            result.total_Bending_energy = total_Bending_energy;
            result.total_Bending_energy_calculated = total_Bending_energy_calculated;
            result.plane_bending_energy = plane_bending_energy;   % NOTE: renamed from plain_Bending_energy
            result.cap_bending_energy   = cap_bending_energy;
            
            result.total_area = total_area;
            result.protein_caps_area = protein_caps_area;
            result.membrane_area = membrane_area;
            result.big_sphere_area = big_sphere_area;
            result.big_sphere_calc_area = big_sphere_calc_area;
            
            result.num_triangles = num_triangles;
            
            result.protein_base_area_A = protein_base_area_A;
            result.protein_base_area_B = protein_base_area_B;
            result.protein_base_area_C = protein_base_area_C;
            result.protein_base_area_D = protein_base_area_D;
            result.protein_base_area_E = protein_base_area_E;
            result.total_protein_base_area = total_protein_base_area;
            
            result.num_proteins = num_proteins;
            result.area_per_protein_anltc = area_per_protein_anltc;
            result.area_per_protein_real  = area_per_protein_real;
            
            result.P_in_model = P_in_model;
            
            result.arc_AB = arc_AB;
            result.arc_AC = arc_AC;
            result.arc_AD = arc_AD;
            result.arc_AE = arc_AE;
            result.arc_BC = arc_BC;
            result.arc_CD = arc_CD;
            result.arc_DE = arc_DE;
            
            result.analytic_total_area = analytic_total_area;
            result.PSI = PSI;
            result.PSI_honey = PSI_val;
            result.PHI = PHI;
            result.RATIO = ratio_value;

            % If the run succeeded, delete the .fe file
            end
    else
        disp('There was an error running Surface Evolver.');
        result.error_message = out_string;  % Log the error message
    end

    % Delete run file only if run was successful
    if status == 0 && exist(fname_run, 'file')
        try delete(fname_run); catch, end
    end

    % delete(fname_run);
end

function fn = getFilename(p), [~,fn,ext] = fileparts(p); fn = [fn '.' ext]; end
function safeDeleteFiles(cellOfPaths)
    for k = 1:numel(cellOfPaths)
        if exist(cellOfPaths{k},'file'), try delete(cellOfPaths{k}); catch,end,end
    end
end
function [idx,val] = fetchNextOrTimeout(futs, open, poll)
idx = []; val = [];
timeoutSec = 60;
try
    [idx,val] = fetchNext(futs, poll);
catch
    idx = find(open & strcmp({futs.State},'finished'),1);
    if ~isempty(idx)
        try val = fetchOutputs(futs(idx)); catch, val = []; end
    end
end

% cancel any still‑running future that exceeded the limit
for k = find(open)
    if seconds(futs(k).RunningDuration) > timeoutSec
        cancel(futs(k)); open(k)=false;
    end
end
end

function S = buildNaNStruct(radP,J1,J2,avgCurv,distVal,nProt)
% Same fields as the real result, but:
%   • geometry / input columns are still filled
%   • everything the worker should have produced is NaN
%
S = struct( ...
  'radius_P1',radP, 'radius_P2',radP, 'radius_P3',radP, 'radius_P4',radP, 'radius_P5',radP, ...
   'Curv_P1_fam', NaN, 'Curv_P2_fam', NaN, ...   
  'Curv_P1',J1, 'Curv_P2',J2, 'Curv_P3',J2, 'Curv_P4',J2, 'Curv_P5',J2, ...
  'Avg_curv',avgCurv, 'distance',distVal, 'membrane_tension',NaN, ...
  'total_energy',NaN, 'total_Bending_energy',NaN, 'total_Bending_energy_calculated',NaN, ...
  'plane_bending_energy',NaN, 'cap_bending_energy',NaN, ...
  'total_area',NaN, 'protein_caps_area',NaN, 'membrane_area',NaN, ...
  'big_sphere_area',NaN, 'big_sphere_calc_area',NaN, ...
  'num_triangles',NaN, ...
  'protein_base_area_A',NaN, 'protein_base_area_B',NaN, 'protein_base_area_C',NaN, ...
  'protein_base_area_D',NaN, 'protein_base_area_E',NaN, ...
  'total_protein_base_area',NaN, ...
  'num_proteins',NaN, 'area_per_protein_anltc',NaN, 'area_per_protein_real',NaN, ...
  'P_in_model',nProt, ...
  'arc_AB',NaN, 'arc_AC',NaN, 'arc_AD',NaN, 'arc_AE',NaN, 'arc_BC',NaN, 'arc_CD',NaN, 'arc_DE',NaN, ...
  'analytic_total_area',NaN, 'PSI',NaN, 'PSI_honey',NaN, 'PHI',NaN,'RATIO',NaN,...
  'error_message','' ...
);

end

function min_Dist = calculate_min_dist(rad_P1, AVG_CURVATURE, NUMBER_OF_PROTEINS)
% CALCULATE_MIN_DIST
% Finds the smallest Dist such that the protein-protein spherical separations
% in the current geometry are all at least 5*rad_P1.
%
% Current geometry:
%   A = protein at (0,0,0)
%   B = protein
%   C = protein
%   D = membrane corner only (not a protein)
%
% Therefore the relevant protein-protein distances are:
%   AB, BC, AC
%
% NUMBER_OF_PROTEINS is kept in the signature for compatibility.

    BIG_SPHERE_RADIUS = 2 / AVG_CURVATURE;
    ZC = -BIG_SPHERE_RADIUS; %#ok<NASGU>

    min_required = 5 * rad_P1;
    tol = 1e-6;
    step = 1;

    Dist = 1;

    while true
        % ============================================================
        % ANCHOR: current SE geometry coordinates
        % ============================================================
        theta = Dist / BIG_SPHERE_RADIUS;
        sinus = sin(theta);
        cosin = cos(theta);
        z_ring = BIG_SPHERE_RADIUS * (cosin - 1);

        % A : central protein
        Ax = 0;
        Ay = 0;
        Az = 0;

        % B : as in SE block
        Bx = BIG_SPHERE_RADIUS * sinus * cos(0);
        By = BIG_SPHERE_RADIUS * sinus * sin(0);
        Bz = z_ring;

        % C : as in SE block
        Cx = -BIG_SPHERE_RADIUS * sinus * cos(pi/3);
        Cy =  BIG_SPHERE_RADIUS * sinus * sin(pi/3);
        Cz = z_ring;

        % D exists geometrically, but is NOT a protein
        % Dx = -BIG_SPHERE_RADIUS * sinus / sqrt(1 + 3 * cosin * cosin);
        % Dy = 0;
        % Dz = BIG_SPHERE_RADIUS * (2 * cosin / sqrt(1 + 3 * cosin * cosin) - 1);

        % ============================================================
        % ANCHOR: protein-protein spherical distances
        % ============================================================
        dAB = calculate_spherical_distance(Ax, Ay, Az, Bx, By, Bz, BIG_SPHERE_RADIUS);
        dBC = calculate_spherical_distance(Bx, By, Bz, Cx, Cy, Cz, BIG_SPHERE_RADIUS);
        dAC = calculate_spherical_distance(Ax, Ay, Az, Cx, Cy, Cz, BIG_SPHERE_RADIUS);

        % Stop once all protein-protein spacings are large enough
        if (dAB >= min_required - tol) && ...
           (dBC >= min_required - tol) && ...
           (dAC >= min_required - tol)

            min_Dist = Dist;
            break
        end

        Dist = Dist + step;
    end

    fprintf('min Dist for current A-B-C protein geometry is: %.3f\n', min_Dist);
end

function distance_value = calculate_spherical_distance( ...
    input_Ax, input_Ay, input_Az, ...
    input_Bx, input_By, input_Bz, ...
    BIG_SPHERE_RADIUS)
% CALCULATE_SPHERICAL_DISTANCE (robust; center-aware; no single letters)
% Geodesic distance along the great circle on a sphere of radius
% BIG_SPHERE_RADIUS whose center is (0,0,ZC).
% Uses atan2(norm(cross(u,v)), dot(u,v)) for numerical stability.

    % ---- locals ----------------------------------------------------
    sphere_radius = BIG_SPHERE_RADIUS;
    ZC = BIG_SPHERE_RADIUS;
    tiny_value = 1e-300;   % protects against 0/0 in extreme degeneracy

    % ---- shift to sphere-centered frame (center at 0,0,ZC) --------
    shifted_Ax = input_Ax;
    shifted_Ay = input_Ay;
    shifted_Az = input_Az - ZC;

    shifted_Bx = input_Bx;
    shifted_By = input_By;
    shifted_Bz = input_Bz - ZC;

    % ---- normalize both points to unit vectors ---------------------
    norm_A_value = sqrt(shifted_Ax.^2 + shifted_Ay.^2 + shifted_Az.^2);
    norm_B_value = sqrt(shifted_Bx.^2 + shifted_By.^2 + shifted_Bz.^2);

    % guard norms (should be ~sphere_radius, but renormalize anyway)
    if norm_A_value < tiny_value, norm_A_value = sphere_radius; end
    if norm_B_value < tiny_value, norm_B_value = sphere_radius; end

    unit_Ax = shifted_Ax ./ norm_A_value;
    unit_Ay = shifted_Ay ./ norm_A_value;
    unit_Az = shifted_Az ./ norm_A_value;

    unit_Bx = shifted_Bx ./ norm_B_value;
    unit_By = shifted_By ./ norm_B_value;
    unit_Bz = shifted_Bz ./ norm_B_value;

    % ---- dot and cross on the unit sphere --------------------------
    dot_unitA_unitB = unit_Ax.*unit_Bx + unit_Ay.*unit_By + unit_Az.*unit_Bz;

    % clamp tiny overshoots to [-1,1]
    dot_unitA_unitB = max(-1, min(1, dot_unitA_unitB));

    cross_x = unit_Ay.*unit_Bz - unit_Az.*unit_By;
    cross_y = unit_Az.*unit_Bx - unit_Ax.*unit_Bz;
    cross_z = unit_Ax.*unit_By - unit_Ay.*unit_Bx;
    cross_norm = sqrt(cross_x.^2 + cross_y.^2 + cross_z.^2);

    % ---- central angle via atan2 (stable for tiny angles) ----------
    central_angle_value = atan2(cross_norm, dot_unitA_unitB);

    % ---- geodesic distance on radius R sphere ----------------------
    distance_value = sphere_radius .* central_angle_value;
end



function range = custom_distance_range(start_val, end_val, protein_number)
    if start_val <= 30
        part1 = start_val:0.25:min(30, end_val);
    else
        part1 = [];
    end

    if end_val > 30 && start_val <= 40
        part2_start = max(start_val, 30.5);
        part2_end = min(end_val, 39.5);
        part2 = part2_start:0.5:part2_end;
    else
        part2 = [];
    end

    if end_val > 40
        part3_start = max(start_val, 40);
        part3 = part3_start:1:end_val;

        if ((end_val > 40) && (protein_number >=11))
        part3_start = max(start_val, 40);
        part3 = part3_start:0.5:end_val;
        end

    else
        part3 = [];
    end

    range = [part1, part2, part3];
end

function angleDeg = calculate_vertex_A_angle( ...
        Ax, Ay, Az, ...
        Bx, By, Bz, ...
        Cx, Cy, Cz)
% CALCULATE_VERTEX_A_ANGLE  Angle ∠BAC on a sphere or in Euclidean space
%
%   angleDeg = CALCULATE_VERTEX_A_ANGLE(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz)
%   takes the Cartesian coordinates of vertices A, B, and C and returns
%   the angle at vertex A (between vectors AB and AC) in degrees.
%
%   Inputs:
%     Ax, Ay, Az – coordinates of vertex A
%     Bx, By, Bz – coordinates of vertex B
%     Cx, Cy, Cz – coordinates of vertex C
%
%   Output:
%     angleDeg   – angle at A, in degrees

    % Vectors AB and AC
    AB = [Bx - Ax, By - Ay, Bz - Az];
    AC = [Cx - Ax, Cy - Ay, Cz - Az];

    % Dot product and magnitudes
    dotProd = dot(AB, AC);
    magAB   = norm(AB);
    magAC   = norm(AC);

    % Cosine of the angle (clip for numerical robustness)
    cosAng = dotProd / (magAB * magAC);
    cosAng = max(-1, min(1, cosAng));

    % Angle in degrees
    angleDeg = rad2deg(acos(cosAng));
end

function area = sector_area(radius, angle_deg)
% SECTOR_AREA  Calculates the area of a sector of a circle
%
%   area = SECTOR_AREA(radius, angle_deg)
%   returns the area of a sector with the given radius and angle (in degrees).
%
%   Inputs:
%     radius     – radius of the circle
%     angle_deg  – central angle in degrees
%
%   Output:
%     area       – sector area

    % Convert angle to radians
    angle_rad = deg2rad(angle_deg);

    % Compute area
    area = 0.5 * radius^2 * angle_rad;
end

function triangle_area_value = calculate_spherical_triangle_area( ...
    Ax, Ay, Az, ...
    Bx, By, Bz, ...
    Cx, Cy, Cz, ...
    BIG_SPHERE_RADIUS)
% CALCULATE_SPHERICAL_TRIANGLE_AREA (Girard/atan2; center-aware; robust)
% Assumes the sphere has radius BIG_SPHERE_RADIUS and center (0,0,ZC).
% Steps:
%   1) Shift vertices by the sphere center.
%   2) Normalize each to the unit sphere.
%   3) Use Girard with atan2: E = 2*atan2(|u·(v×w)|, 1 + u·v + v·w + w·u).
%   4) Area = R^2 * E.
%
% Inputs are scalars (coordinates of A, B, C), plus BIG_SPHERE_RADIUS and ZC.

    % ---- constants / locals ---------------------------------------
    sphere_radius = BIG_SPHERE_RADIUS;
    ZC = BIG_SPHERE_RADIUS;
    tiny_value = 1e-300;   % numerical floor (avoids 0/0 in extreme degeneracy)

    % ---- shift vertices by sphere center (0,0,ZC) ------------------
    shifted_Ax = Ax;
    shifted_Ay = Ay;
    shifted_Az = Az - ZC;

    shifted_Bx = Bx;
    shifted_By = By;
    shifted_Bz = Bz - ZC;

    shifted_Cx = Cx;
    shifted_Cy = Cy;
    shifted_Cz = Cz - ZC;

    % ---- normalize to unit sphere ---------------------------------
    norm_A_value = sqrt(shifted_Ax.^2 + shifted_Ay.^2 + shifted_Az.^2);
    norm_B_value = sqrt(shifted_Bx.^2 + shifted_By.^2 + shifted_Bz.^2);
    norm_C_value = sqrt(shifted_Cx.^2 + shifted_Cy.^2 + shifted_Cz.^2);

    % guard norms (should be ~sphere_radius)
    if norm_A_value < tiny_value, norm_A_value = sphere_radius; end
    if norm_B_value < tiny_value, norm_B_value = sphere_radius; end
    if norm_C_value < tiny_value, norm_C_value = sphere_radius; end

    unit_Ax = shifted_Ax ./ norm_A_value;
    unit_Ay = shifted_Ay ./ norm_A_value;
    unit_Az = shifted_Az ./ norm_A_value;

    unit_Bx = shifted_Bx ./ norm_B_value;
    unit_By = shifted_By ./ norm_B_value;
    unit_Bz = shifted_Bz ./ norm_B_value;

    unit_Cx = shifted_Cx ./ norm_C_value;
    unit_Cy = shifted_Cy ./ norm_C_value;
    unit_Cz = shifted_Cz ./ norm_C_value;

    % ---- dot products on unit sphere -------------------------------
    dot_A_B = unit_Ax.*unit_Bx + unit_Ay.*unit_By + unit_Az.*unit_Bz;
    dot_B_C = unit_Bx.*unit_Cx + unit_By.*unit_Cy + unit_Bz.*unit_Cz;
    dot_C_A = unit_Cx.*unit_Ax + unit_Cy.*unit_Ay + unit_Cz.*unit_Az;

    % clamp to [-1,1] (safety against tiny overshoots)
    dot_A_B = max(-1, min(1, dot_A_B));
    dot_B_C = max(-1, min(1, dot_B_C));
    dot_C_A = max(-1, min(1, dot_C_A));

    % ---- cross product and triple scalar product -------------------
    cross_BC_x = unit_By.*unit_Cz - unit_Bz.*unit_Cy;
    cross_BC_y = unit_Bz.*unit_Cx - unit_Bx.*unit_Cz;
    cross_BC_z = unit_Bx.*unit_Cy - unit_By.*unit_Cx;

    triple_A_BC = unit_Ax.*cross_BC_x + unit_Ay.*cross_BC_y + unit_Az.*cross_BC_z;

    % ---- Girard (robust atan2) on unit sphere ----------------------
    denominator_value = 1 + dot_A_B + dot_B_C + dot_C_A;
    spherical_excess_value = 2 .* atan2(abs(triple_A_BC), denominator_value);

    % ---- final area on radius R sphere ------------------------------
    triangle_area_value = sphere_radius.^2 .* spherical_excess_value;
end




function [PSI, PHI, membrane_area] = calculate_PSI(Rad_P, vertices, RATIO_VAR, BIG_SPHERE_RADIUS_var)
% CALCULATE_PSI for the current geometry:
% - A, B, C are proteins
% - D is only a membrane corner, not a protein
%
% Sector definitions:
% - Sector at A is the sum of two wedges around A:
%       (A_B, A_C) + (A_C, A_D)
% - Sector at B is (B_A, B_C)
% - Sector at C is (C_B, C_D)
% - D contributes no protein base area
%
% Membrane patch:
% - area(BCD)
%
% Returns:
%   PSI           : currently assigned from RATIO_VAR exactly as in your old code
%   PHI           : total protein base coverage / membrane patch area
%   membrane_area : area(BCD)

% --------- sanity check: required fields ----------
req = {'A','B','C','D', ...
       'A_B','A_C','A_D', ...
       'B_A','B_C', ...
       'C_B','C_D'};
for k = 1:numel(req)
    if ~isfield(vertices, req{k})
        error('calculate_PSI:MissingField', 'vertices.%s is missing.', req{k});
    end
end

% --------- helper for an angle at center P using two periphery points U,V ----------
ang = @(P, U, V) calculate_vertex_A_angle( ...
    P(1), P(2), P(3), ...
    U(1), U(2), U(3), ...
    V(1), V(2), V(3));

% --------- sector angles (radians) ----------
angle_A = ang(vertices.A, vertices.A_B, vertices.A_C) ...
        + ang(vertices.A, vertices.A_C, vertices.A_D);

angle_B = ang(vertices.B, vertices.B_A, vertices.B_C);
angle_C = ang(vertices.C, vertices.C_B, vertices.C_D);

% --------- base areas (sectors on the small bases) ----------
area_A = sector_area(Rad_P, angle_A);
area_B = sector_area(Rad_P, angle_B);
area_C = sector_area(Rad_P, angle_C);

totalBase = area_A + area_B + area_C;

% --------- PSI ----------
if RATIO_VAR == 1
    PSI = 1;
elseif RATIO_VAR == 11
    PSI = 0.5;
elseif RATIO_VAR == 13
    PSI = 0.25;
elseif RATIO_VAR == 15
    PSI = 0.17;
else
    PSI = NaN;
end

% --------- membrane patch area = triangle BCD ----------
membrane_area = calculate_spherical_triangle_area( ...
    vertices.B(1), vertices.B(2), vertices.B(3), ...
    vertices.C(1), vertices.C(2), vertices.C(3), ...
    vertices.D(1), vertices.D(2), vertices.D(3), ...
    BIG_SPHERE_RADIUS_var);

% --------- PHI: total protein base coverage over the membrane patch ----------
PHI = totalBase ./ membrane_area;

% optional diagnostics
% fprintf('membrane area is = %.8f\n', membrane_area);
% disp(PSI), disp(PHI)
end

function membrane_area = calculate_AREA(vertices, BIG_SPHERE_RADIUS_var)
% CALCULATE_AREA
% Returns only the spherical area of triangle BCD.
%
% IMPORTANT:
%   vertices.B, vertices.C, vertices.D must be the actual mesh vertices.

% --------- sanity check: required fields ----------
req = {'B','C','D'};
for k = 1:numel(req)
    if ~isfield(vertices, req{k})
        error('calculate_AREA:MissingField', 'vertices.%s is missing.', req{k});
    end
end

% --------- membrane patch area = triangle BCD ----------
membrane_area = calculate_spherical_triangle_area( ...
    vertices.B(1), vertices.B(2), vertices.B(3), ...
    vertices.C(1), vertices.C(2), vertices.C(3), ...
    vertices.D(1), vertices.D(2), vertices.D(3), ...
    BIG_SPHERE_RADIUS_var);

% fprintf('membrane area = %.8f\n', membrane_area);
end


function vertices = get_model_vertices(Dist_var, P_var, AVG_CURVATURE_var, NUMBER_OF_PROTEINS_var)
% GET_MODEL_VERTICES
% Mirrors the CURRENT Surface Evolver parameter block exactly:
%
% - A is the center cap at (0,0,0)
% - B is on the ring at azimuth 0
% - C is on the ring at the SE-defined position:
%       Cx = -R*sin(theta)*cos(pi/3)
%       Cy =  R*sin(theta)*sin(pi/3)
%       Cz = z_ring
% - D is the special lower point defined by:
%       Dx = -R*sin(theta)/sqrt(1+3*cos(theta)^2)
%       Dy = 0
%       Dz = R*(2*cos(theta)/sqrt(1+3*cos(theta)^2) - 1)
%
% The function returns:
% - cap points A,B,C,D
% - bases BaseA..BaseD
% - cap tops A_top..D_top
% - cap centers CtrA..CtrD
% - triangle centers TABC, TACD
% - periphery points exactly as in the SE block
%
% NOTE:
% In the SE "vertices" section, the actual mesh vertex for D is D = [Dx Dy Dz],
% not D_top. This function keeps both:
%   v.D     = [Dx Dy Dz]
%   v.D_top = cap-top point from the parameter block

% ============================================================
% ANCHOR: inputs -> model params
% ============================================================
Dist               = Dist_var;
rad_P1             = P_var;
rad_P2             = P_var;
rad_P3             = P_var;
rad_P4             = P_var;

NUMBER_OF_PROTEINS = NUMBER_OF_PROTEINS_var;

spin_angle_deg     = 360 / 6;
spin_angle         = deg2rad(spin_angle_deg); %#ok<NASGU>

AVG_CURVATURE      = AVG_CURVATURE_var;
BIG_SPHERE_RADIUS  = 2 / AVG_CURVATURE;
ZC                 = -BIG_SPHERE_RADIUS;

% ============================================================
% ANCHOR: protein curvature block
% ============================================================
Average_curvature_P_general = 0.01;

Average_curvature_P1 = Average_curvature_P_general;
Average_curvature_P2 = Average_curvature_P_general;
Average_curvature_P3 = Average_curvature_P_general;
Average_curvature_P4 = Average_curvature_P_general;

Rcap_P1 = 2 / Average_curvature_P1;
Rcap_P2 = 2 / Average_curvature_P2;
Rcap_P3 = 2 / Average_curvature_P3;
Rcap_P4 = 2 / Average_curvature_P4;

Hcap_P1 = Rcap_P1 - sqrt(Rcap_P1^2 - rad_P1^2);
Hcap_P2 = Rcap_P2 - sqrt(Rcap_P2^2 - rad_P2^2);
Hcap_P3 = Rcap_P3 - sqrt(Rcap_P3^2 - rad_P3^2);
Hcap_P4 = Rcap_P4 - sqrt(Rcap_P4^2 - rad_P4^2);

% ============================================================
% ANCHOR: big-sphere shorthand
% ============================================================
R     = BIG_SPHERE_RADIUS;
theta = Dist / BIG_SPHERE_RADIUS;
sinus = sin(theta);
cosin = cos(theta);
z_ring = BIG_SPHERE_RADIUS * (cosin - 1);

sphere_center = [0, 0, ZC];

% ============================================================
% ANCHOR: cap positions A..D (EXACTLY as in current SE block)
% ============================================================
A = [0, 0, 0];

B = [ ...
    BIG_SPHERE_RADIUS * sinus * cos(0), ...
    BIG_SPHERE_RADIUS * sinus * sin(0), ...
    z_ring];

C = [ ...
   -BIG_SPHERE_RADIUS * sinus * cos(pi/3), ...
    BIG_SPHERE_RADIUS * sinus * sin(pi/3), ...
    z_ring];

D = [ ...
   -BIG_SPHERE_RADIUS * sinus / sqrt(1 + 3 * cosin * cosin), ...
    0, ...
    BIG_SPHERE_RADIUS * (2 * cosin / sqrt(1 + 3 * cosin * cosin) - 1)];

% ============================================================
% ANCHOR: bases (contact circles)
% ============================================================
PlaneOffset = @(rad) BIG_SPHERE_RADIUS * sqrt(1 - (rad / BIG_SPHERE_RADIUS)^2);
BaseScale   = @(rad) PlaneOffset(rad) / BIG_SPHERE_RADIUS;

BaseA = [ ...
    A(1), ...
    A(2), ...
    ZC + (A(3) - ZC) * BaseScale(rad_P1)];

BaseB = [ ...
    B(1) * BaseScale(rad_P2), ...
    B(2) * BaseScale(rad_P2), ...
    ZC + (B(3) - ZC) * BaseScale(rad_P2)];

BaseC = [ ...
    C(1) * BaseScale(rad_P3), ...
    C(2) * BaseScale(rad_P3), ...
    ZC + (C(3) - ZC) * BaseScale(rad_P3)];

BaseD = [ ...
    D(1) * BaseScale(rad_P4), ...
    D(2) * BaseScale(rad_P4), ...
    ZC + (D(3) - ZC) * BaseScale(rad_P4)];

% ============================================================
% ANCHOR: unit directions for cap axes
% ============================================================
unit_dir = @(Base) (Base - sphere_center) / norm(Base - sphere_center);

dirA = unit_dir(BaseA);
dirB = unit_dir(BaseB);
dirC = unit_dir(BaseC);
dirD = unit_dir(BaseD);

% ============================================================
% ANCHOR: cap tops
% ============================================================
Atop = BaseA + dirA * Hcap_P1;
Btop = BaseB + dirB * Hcap_P2;
Ctop = BaseC + dirC * Hcap_P3;
Dtop = BaseD + dirD * Hcap_P4;

% ============================================================
% ANCHOR: small-sphere centres
% ============================================================
CtrA = BaseA - dirA * (Rcap_P1 - Hcap_P1);
CtrB = BaseB - dirB * (Rcap_P2 - Hcap_P2);
CtrC = BaseC - dirC * (Rcap_P3 - Hcap_P3);
CtrD = BaseD - dirD * (Rcap_P4 - Hcap_P4);

% ============================================================
% ANCHOR: renormalization to big sphere
% Matches your SE formula:
%   scale = R / sqrt(x^2 + y^2 + (z-ZC)^2)
%   P_out = scale * P
% ============================================================
renorm = @(P) (R / sqrt(P(1)^2 + P(2)^2 + (P(3) - ZC)^2)) * P;

% ============================================================
% ANCHOR: triangle centers
% ============================================================
tABC_1 = (A + B + C) / 3;
tACD_1 = (A + C + D) / 3;

TABC = renorm(tABC_1);
TACD = renorm(tACD_1);

% ============================================================
% ANCHOR: periphery builder
% ============================================================
periph = @(BaseP, Target, rad_P) renorm( ...
    BaseP + (rad_P / norm(Target - BaseP)) * (Target - BaseP) );

% ============================================================
% ANCHOR: A periphery
% ============================================================
A_B    = periph(BaseA, B,    rad_P1);
A_C    = periph(BaseA, C,    rad_P1);
A_TABC = periph(BaseA, TABC, rad_P1);
A_D    = periph(BaseA, D,    rad_P1);
A_TACD = periph(BaseA, TACD, rad_P1);

% ============================================================
% ANCHOR: B periphery
% ============================================================
B_A    = periph(BaseB, A,    rad_P2);
B_C    = periph(BaseB, C,    rad_P2);
B_TABC = periph(BaseB, TABC, rad_P2);

% ============================================================
% ANCHOR: C periphery
% ============================================================
C_B    = periph(BaseC, B,    rad_P3);
C_A    = periph(BaseC, A,    rad_P3);
C_TABC = periph(BaseC, TABC, rad_P3);
C_D    = periph(BaseC, D,    rad_P3);
C_TACD = periph(BaseC, TACD, rad_P3);

% ============================================================
% ANCHOR: D periphery
% ============================================================
D_C    = periph(BaseD, C,    rad_P4);
D_A    = periph(BaseD, A,    rad_P4);
D_TACD = periph(BaseD, TACD, rad_P4);

% ============================================================
% ANCHOR: pack output
% ============================================================
v = struct();

% Big-sphere cap positions / mesh points
v.A = A;
v.B = B;
v.C = C;
v.D = D;

% Bases
v.BaseA = BaseA;
v.BaseB = BaseB;
v.BaseC = BaseC;
v.BaseD = BaseD;

% Cap tops
v.A_top = Atop;
v.B_top = Btop;
v.C_top = Ctop;
v.D_top = Dtop;

% Small-sphere centres
v.CtrA = CtrA;
v.CtrB = CtrB;
v.CtrC = CtrC;
v.CtrD = CtrD;

% Triangle centers
v.TABC = TABC;
v.TACD = TACD;

% A periphery
v.A_B    = A_B;
v.A_C    = A_C;
v.A_TABC = A_TABC;
v.A_D    = A_D;
v.A_TACD = A_TACD;

% B periphery
v.B_A    = B_A;
v.B_C    = B_C;
v.B_TABC = B_TABC;

% C periphery
v.C_B    = C_B;
v.C_A    = C_A;
v.C_TABC = C_TABC;
v.C_D    = C_D;
v.C_TACD = C_TACD;

% D periphery
v.D_C    = D_C;
v.D_A    = D_A;
v.D_TACD = D_TACD;

% Bookkeeping
v.R = R;
v.ZC = ZC;
v.theta = theta;
v.sinus = sinus;
v.cosin = cosin;
v.z_ring = z_ring;
v.spin_angle_deg = spin_angle_deg;
v.NUMBER_OF_PROTEINS = NUMBER_OF_PROTEINS;

% Protein parameters
v.rad_P1 = rad_P1;
v.rad_P2 = rad_P2;
v.rad_P3 = rad_P3;
v.rad_P4 = rad_P4;

v.Rcap_P1 = Rcap_P1;
v.Rcap_P2 = Rcap_P2;
v.Rcap_P3 = Rcap_P3;
v.Rcap_P4 = Rcap_P4;

v.Hcap_P1 = Hcap_P1;
v.Hcap_P2 = Hcap_P2;
v.Hcap_P3 = Hcap_P3;
v.Hcap_P4 = Hcap_P4;

vertices = v;
end



% ---- Helper: Slerp midpoint ----
function [x, y, z] = slerp_midpoint(Bx, By, Bz, Cx, Cy, Cz, R, ZC)
    B = [Bx; By; Bz - ZC]; C = [Cx; Cy; Cz - ZC];
    Bu = B / norm(B); Cu = C / norm(C);
    omega = acos(dot(Bu, Cu));
    sin_omega = sin(omega);
    scaleB = sin(0.5 * omega) / sin_omega;
    scaleC = scaleB;
    mid = R * (scaleB * Bu + scaleC * Cu);
    x = mid(1);
    y = mid(2);
    z = mid(3) + ZC;
end

% ---- Helper: radial point on sphere ----
function xyz = spherical_radial_point(BaseX, BaseY, BaseZ, Tx, Ty, Tz, r, R, ZC)
    Nx = r / sqrt((Tx - BaseX)^2 + (Ty - BaseY)^2 + (Tz - BaseZ)^2);
    Px = BaseX + (Tx - BaseX) * Nx;
    Py = BaseY + (Ty - BaseY) * Nx;
    Pz = BaseZ + (Tz - BaseZ) * Nx;
    scale = R / sqrt(Px^2 + Py^2 + (Pz - ZC)^2);
    xyz = [Px * scale, Py * scale, Pz * scale];
end

% ---- Helper: unit direction from base to center ----
function dir = unit_direction(x, y, z, ZC)
    dz = z - ZC;
    norm_val = sqrt(x^2 + y^2 + dz^2);
    dir = [x, y, dz] / norm_val;
end

function stats = run_model(Min_Dist_var, Max_Dist_var, ...
                           Prot_rad_var, AVG_CURVATURE_var, ...
                           RATIO_VAR, num_PHI_values, make_plots)
% RUN_MODEL
%   Builds a Dist grid with exactly num_PHI_values samples per N,
%   biased to be DENSE at LOW area and Sparser at HIGH area.
%   Also logs PSI/PHI ranges and produces diagnostics plots (optional).
%
% REQUIREMENT:
%   [PSI, PHI, membrane_area] = calculate_PSI(Rad_P, vertices, N, BIG_SPHERE_RADIUS)
%
% INPUTS:
%   Min_Dist_var, Max_Dist_var  : requested Dist bounds (nm)
%   Prot_rad_var                : protein base radius (nm)
%   AVG_CURVATURE_var           : big-sphere average curvature (nm^-1)
%   NUMBER_OF_PROTEINS          : vector of N values (e.g., [0 1 2 3 4 5])
%   num_PHI_values              : target number of Dist samples per N
%   make_plots                  : (optional) true/false (default: true)
%
% OUTPUT:
%   stats(k) with fields:
%     .Protein_rad, .PROTEINS
%     .PSI_min/.PSI_max, .PHI_min/.PHI_max
%     .Dist_grid  (1×num_PHI_values)

if nargin < 7 || isempty(make_plots), make_plots = true; end

% ----------------------------- constants / knobs ------------------------
BIG_SPHERE_RADIUS = 2 / AVG_CURVATURE_var;
numTarget         = num_PHI_values;   % samples per N
denseSweepPts     = 401;              % resolution for area(Dist) sweep
gamma             = 1;                % weighting power: 1 (mild) … 2 (stronger)

% ----------------------------- preallocate output -----------------------
stats = struct( ...
    'Protein_rad', Prot_rad_var, ...
    'PROTEINS',    [], ...
    'PSI_min',     [], 'PSI_max', [], ...
    'PHI_min',     [], 'PHI_max', [], ...
    'Dist_grid',   [] );

% --------------------------------- main loop ----------------------------
for k = 1:numel(RATIO_VAR)
    N = RATIO_VAR(k);

    % --- geometric minimum distance ---
    geomMin = calculate_min_dist(Prot_rad_var, AVG_CURVATURE_var, N);
    Dist_lo = max(Min_Dist_var, geomMin);
    Dist_hi = max(Dist_lo, Max_Dist_var);   % ensure hi ≥ lo

    % --- dense sweep: Dist → membrane area(Dist) ---
    Dist_sweep = linspace(Dist_lo, Dist_hi, denseSweepPts);
    A_mem      = nan(size(Dist_sweep));

    for j = 1:numel(Dist_sweep)
        Dist = Dist_sweep(j);
        v = get_model_vertices(Dist, Prot_rad_var, AVG_CURVATURE_var, N);
        [A_mem(j)] = calculate_AREA(v, BIG_SPHERE_RADIUS);
    end

    % keep valid, positive areas
    valid = isfinite(A_mem) & (A_mem > 0);
    Ds_all = Dist_sweep(valid);
    As_all = A_mem(valid);

    % ------------------ build Dist grid (weighted-arc) -------------------
    if numel(Ds_all) < 2
        warning('N=%d: insufficient valid area samples; using linear Dist grid.', N);
        Dist_grid      = linspace(Dist_lo, Dist_hi, numTarget);
        sqrtA_on_grid  = nan(size(Dist_grid));   % plotting fallback
    else
        % Weight: larger when area is small → denser sampling at low area
        w = 1 ./ (As_all .^ gamma);

        % Cumulative “arc” in weighted space
        C = cumtrapz(Ds_all, w);

        % Ensure strictly increasing for interpolation
        [Cu, ia] = unique(C, 'stable');
        Du       = Ds_all(ia);

        if numel(Cu) < 2
            warning('N=%d: degenerate cumulative; using linear Dist grid.', N);
            Dist_grid     = linspace(Dist_lo, Dist_hi, numTarget);
            sqrtA_on_grid = nan(size(Dist_grid));
        else
            % Equal steps in cumulative → map back to Dist
            C_targets = linspace(Cu(1), Cu(end), numTarget);
            Dist_grid = interp1(Cu, Du, C_targets, 'pchip');

            % numeric safety / clamp to bounds
            Dist_grid = real(Dist_grid(:).');                % row vector
            Dist_grid = max(min(Dist_grid, Dist_hi), Dist_lo);

            % compute sqrt(area) at selected Dist (for overlay plot)
            sqrtA_on_grid = nan(size(Dist_grid));
            for j = 1:numel(Dist_grid)
                Dist = Dist_grid(j);
                v = get_model_vertices(Dist, Prot_rad_var, AVG_CURVATURE_var, N);
                [A_here] = calculate_AREA(v, BIG_SPHERE_RADIUS);
                sqrtA_on_grid(j) = sqrt(A_here);
            end
        end
    end

    % --- evaluate PSI/PHI on the final grid (for stats only) ---
    PSI_vals = zeros(size(Dist_grid));
    PHI_vals = zeros(size(Dist_grid));
    for j = 1:numel(Dist_grid)
        Dist = Dist_grid(j);
        v = get_model_vertices(Dist, Prot_rad_var, AVG_CURVATURE_var, N);
        [PSI_vals(j), PHI_vals(j)] = ...
            calculate_PSI(Prot_rad_var, v, N, BIG_SPHERE_RADIUS);
    end

    % --- stash results ---
    stats(k).Protein_rad = Prot_rad_var;
    stats(k).PROTEINS    = N;
    stats(k).PSI_min     = min(PSI_vals);
    stats(k).PSI_max     = max(PSI_vals);
    stats(k).PHI_min     = min(PHI_vals);
    stats(k).PHI_max     = max(PHI_vals);
    stats(k).Dist_grid   = Dist_grid;   % exactly numTarget points

    fprintf(['N=%d  Dist[%.2f–%.2f]  (%d pts, weighted low→high area, \\gamma=%g)  ', ...
             '| PSI[min=%.4f max=%.4f]  PHI[min=%.4f max=%.4f]\n'], ...
            N, Dist_grid(1), Dist_grid(end), numTarget, gamma, ...
            stats(k).PSI_min, stats(k).PSI_max, ...
            stats(k).PHI_min, stats(k).PHI_max);

    % --------------------------- plots ---------------------------
    if make_plots && numel(Ds_all) >= 2
        % Plot 1: sqrt(area) vs Dist, with selected points overlaid
        figure('Name', sprintf('sqrt(area) vs Dist (N=%d)', N));
        plot(Ds_all, sqrt(As_all), 'LineWidth', 1.5); hold on;
        if all(isfinite(sqrtA_on_grid))
            plot(Dist_grid, sqrtA_on_grid, 'o', 'MarkerSize', 6, 'LineWidth', 1.2);
            legend({'\surd area(Dist) (dense sweep)','selected points'}, 'Location','best');
        else
            legend({'\surd area(Dist) (dense sweep)'}, 'Location','best');
        end
        xlabel('Dist'); ylabel('sqrt(membrane\_area)');
        title(sprintf('N = %d: weighted sampling (\\gamma = %g) → %d points', N, gamma, numTarget));
        grid on; box on;

        % Plot 2: ΔDist spacing — should start small and grow
        dDist = diff(Dist_grid);
        figure('Name', sprintf('ΔDist spacing (N=%d)', N));
        stem(1:numel(dDist), dDist, 'filled');
        xlabel('interval index'); ylabel('\Delta Dist');
        title(sprintf('Gaps between consecutive Dist (N=%d, \\gamma = %g)', N, gamma));
        grid on; box on;

        if ~isempty(dDist)
            fprintf('N=%d: ΔDist — median=%.4f, min=%.4f, max=%.4f\n', ...
                N, median(dDist), min(dDist), max(dDist));
        end
    end
end
end




function results = analyzeProteinCurvature(stats, protein_number, Ja, Jb)
% ANALYZEPROTEINCURVATURE  –  Compute Jint and a 20‑point J range grid
%                             using only the information inside the
%                             `stats` struct produced by run_model.
%
% INPUTS
%   stats          – struct array with fields
%                       .PROTEINS
%                       .PSI_min / PSI_max
%                       .PHI_min / PHI_max
%                       .Dist_grid           (1×N)  ⇐ distance grid
%   protein_number – scalar (e.g. 6, 9, 12 …)
%   Ja, Jb         – curvature values for P1 and P2/3
%
% OUTPUT
%   results struct with fields
%       .PHI_values      1×N     (descending from φmax → φmin)
%       .Jint_values     1×N     (centre curvature for each φ)
%       .J_range_values  N×20    (linear samples between parabola roots)
%       .PSI             scalar  (mean PSI for this protein count)

% ----------------------------------------------------------------------
% 1) Locate the entry for this protein count
% ----------------------------------------------------------------------
idx = find([stats.PROTEINS] == protein_number, 1);
if isempty(idx)
    error('Protein count %d not found in stats struct.', protein_number);
end

PSI      = 0.5 * (stats(idx).PSI_min + stats(idx).PSI_max);
phiRange = [stats(idx).PHI_min, stats(idx).PHI_max];    % [φmin φmax]
distGrid = stats(idx).Dist_grid;                        % 1×N
N        = numel(distGrid);                             % # distance points

% ----------------------------------------------------------------------
% 2) Build the PHI vector (descending – matches original convention)
% ----------------------------------------------------------------------
PHI_values = linspace(phiRange(2), phiRange(1), N);

% ----------------------------------------------------------------------
% 3) Compute Jint and J‑range for each φ
% ----------------------------------------------------------------------
split          = 20;                     % points per row
Jint_values    = zeros(1, N);
J_range_values = zeros(N, split);

PSI = 1-PSI;

for i = 1:N
    phi  = PHI_values(i);
    Jint = phi * (PSI * Ja + (1-PSI) * Jb);
    Jint_values(i) = Jint;

    % Diagnostics (kept from earlier version)
    fprintf('Jint = %.8f\n', Jint);
    fprintf('phi  = %.8f\n', phi);
    fprintf('PSI  = %.8f\n', PSI);

    % Parabola roots and 20 linear samples between them
    [X1, X2] = solve_parabola_vertex_form(Ja, Jb, PSI, Jint);
    J_range_values(i, :) = linspace(X1, X2, split);
end

% ----------------------------------------------------------------------
% 4) Return results
% ----------------------------------------------------------------------
results = struct( ...
    'PHI_values',     PHI_values, ...
    'Jint_values',    Jint_values, ...
    'J_range_values', J_range_values, ...
    'PSI',            PSI );
end

function [left_edge, right_edge] = solve_parabola_vertex_form(Ja, Jb, PSI, Jint)
%SOLVE_PARABOLA_VERTEX_FORM  Roots of  y = width·(x–Jint)^2 + Fmin
%                            with the horizontal line  y = Fmin·multFactor.
%
% INPUTS
%   width       positive scalar  (parabola curvature, “a”)
%   Fmin        vertex height    (≥ 0 in typical use)
%   Jint        vertex x-position
%   multFactor  ≥ 1  ➜  target height = Fmin*multFactor
%   doPlot      (optional, default=false)  draw a figure if true
%
% OUTPUTS
%   left_edge   x-coordinate of left intersection  (< Jint)
%   right_edge  x-coordinate of right intersection (> Jint)
%
% Y. Navot · 2025-07
% -------------------------------------------------------------------------

high_Jcap = max(Ja, Jb);
low_Jcap = min(Ja, Jb);
Delta_Jcap = high_Jcap - low_Jcap;

    % -- algebra ----------------------------------------------------------
    delta = Jint;
    left_edge  = Jint - delta*0.9999;

    if (Ja == low_Jcap) && (Ja ~= Jb) && (PSI ~= 0.5) && (Delta_Jcap >= 0.11)
        delta = 3*delta;
    end

    right_edge = Jint + delta;    
    fprintf('Jint = %.8f\n', Jint);
    fprintf('Roots     : left = %.10g ,  right = %.10g\n', ...
            left_edge, right_edge);
end

% function run_model (Min_Dist_var, Max_Dist_var, Prot_rad_var, AVG_CURVATURE_var, NUMBER_OF_PROTEINS)
% BIG_SPHERE_RADIUS = 2 / AVG_CURVATURE_var;
% 
% min_Dist = calculate_min_dist(Prot_rad_var, AVG_CURVATURE_var, NUMBER_OF_PROTEINS);
% 
% fprintf('min Dist for %d proteins is: %.3f\n', NUMBER_OF_PROTEINS, min_Dist);
% if (Min_Dist_var < min_Dist)
%     Min_Dist_var = min_Dist;
% end
% 
% vertices = get_model_vertices(Min_Dist_var, Prot_rad_var, AVG_CURVATURE_var, NUMBER_OF_PROTEINS);
% 
% [PSI, PHI] = calculate_PSI (Prot_rad_var, vertices, BIG_SPHERE_RADIUS);
% 
% end