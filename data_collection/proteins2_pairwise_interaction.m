close all;
clear all;
clc;

%% this script is for running with larger PHI of all proteins


% Jcap grids as in your current script
% Jcap1 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];
% Jcap2 = [0.001, 0.016, 0.046, 0.076, 0.105, 0.135];

% ----- geometry / model parameters (UNCHANGED from your script) ------
% work PC
% template_dir = 'D:\yiftach_OneDrive\OneDrive - Technion\work_folder\clean scripts\hexagon\data collection\template for surface evolver\actual run templates\';  % template folder
% output_dir   = 'D:\yiftach_OneDrive\OneDrive - Technion\work_folder\my_output_directory_parallel\pairwise_interaction\V1\';        % results root
% kaiser
template_dir = 'E:\Yiftach\OneDrive - Technion\work_folder\clean scripts\hexagon\data collection\template for surface evolver\actual run templates\';  % template folder
output_dir   = 'E:\Yiftach\OneDrive - Technion\kaiser\my_output_directory_parallel\pairwise_interactions\V4\';        % results root

% height_vec = [0, 2, 4, 6, 8, 10];
% dist_vec = [25, 30, 35, 40, 45, 50];
height_vec = [0, 0.7, 1.4, 2];
dist_vec = [25, 30, 35, 40, 45];
Jcap1 = 0.135;
Jcap2 = 0.001;
max_tilt_angle = 9;

% result = run_single_evolver_iteration(0, 40, Jcap1, Jcap2, 25, template_dir, output_dir);

% result = run_single_evolver_iteration(Height_val, Distance, Jcap1, Jcap2, tilt_val, template_dir, output_dir);

run_flat_tilt_scan(Jcap1, Jcap2, dist_vec, height_vec, max_tilt_angle, template_dir, output_dir)

function run_flat_tilt_scan(Jcap1, Jcap2, dist_vec, height_vec, max_tilt_angle, template_dir, output_dir)
% RUN_FLAT_TILT_SCAN(dist_vec, height_vec, max_tilt_angle)
% ---------------------------------------------------------
% dist_vec      : vector of distances [N]
% height_vec    : vector of heights   [M]
% max_tilt_angle: maximal tilt (same units as used in the .fe template)
%
% Inside, all your old machinery (run_model, stats, Jcap loops, ratios,
% parallel watchdogs, etc.) is still used. The only structural change is
% that the inner Evolver calls now loop like:
%
%   for dist(i)
%     for height(j)
%       parfor tilt(k)
%         run_single_evolver_iteration(...)
%       end
%     end
%   end
%
% Save this file as run_flat_tilt_scan.m and call it like:
%   run_flat_tilt_scan([25 30 35], [0 2 4], 30);

    close all;
    clc;

    % ----- defaults if arguments are omitted -----------------------------
    if nargin < 1 || isempty(dist_vec)
        % <<< EDIT this to your preferred default distances >>>
        dist_vec = 25:5:45;     % example: [25 30 35 40 45]
    end
    if nargin < 2 || isempty(height_vec)
        % <<< EDIT this to your preferred default heights >>>
        height_vec = 0;         % example: [0 2 4]
    end
    if nargin < 3 || isempty(max_tilt_angle)
        % <<< EDIT this to your preferred max tilt (deg or rad, as in template) >>>
        max_tilt_angle = 30;
    end

    Prot_rad_var      = 5;
    make_plots        = false;
    sim_only          = false;


    % ----- launch all Jcap / ratio / protein runs ------------------------
    check_curvature_combinations( ...
        Prot_rad_var, ...
        make_plots, ...
        Jcap1, Jcap2, output_dir, template_dir, ...
        sim_only, ...
        dist_vec, height_vec, max_tilt_angle);

    disp('All iterations complete – MAT tables written.');
end


function check_curvature_combinations( ...
    Prot_rad_var, make_plots, ...
    Jcap1, Jcap2, results_folder, template_dir, ...
    sim_only, dist_vec, height_vec, max_tilt_angle)

% CHECK_CURVATURE_COMBINATIONS  – no RATIO, no NUMBER_OF_PROTEINS.
%
% Loops over all (J1,J2) in Jcap1×Jcap2 and, for each pair:
%   • optionally skips if sim_only==true and J1 != J2
%   • creates a pair-root folder
%   • calls activate_and_monitor_parallel(...) which will itself
%     skip existing height-level folders.
%
% Folder layout PER PAIR:
%   results_folder/
%       results_curvatures_[J1, J2]/      <- pair_root
%           results_Dist_xxx/
%               results_Height_yyy/       <- height-level skip is here
%                   results_*.mat
%
% This function does NOT scan over different protein counts anymore.

    %#ok<*NASGU>  % make_plots not used here but kept for signature compatibility

    if nargin < 7 || isempty(sim_only)
        sim_only = false;
    end
    validateattributes(sim_only, {'logical','numeric'}, {'scalar'}, mfilename, 'sim_only');
    sim_only = logical(sim_only);

    if isempty(Jcap1) || isempty(Jcap2)
        error('Jcap1 and Jcap2 must not be empty');
    end

    tol = 1e-12;

    for j = 1:numel(Jcap1)
        for k = 1:numel(Jcap2)
            J1 = Jcap1(j);
            J2 = Jcap2(k);

            J1_str = sprintf('%.3f', J1);
            J2_str = sprintf('%.3f', J2);

            % -------- optional filter for "similar-only" mode --------
            if sim_only && abs(J1 - J2) > tol
                fprintf('• sim_only=TRUE → skip unequal pair [%s,%s].\n', J1_str, J2_str);
                continue;
            end

            % One folder per curvature pair – now prefixed with "results_"
            pair_folder_name = sprintf('results_curvatures_[%s, %s]', J1_str, J2_str);
            pair_root        = fullfile(results_folder, pair_folder_name);

            % Ensure pair folder exists (no skipping here anymore)
            if ~isfolder(pair_root)
                try
                    mkdir(pair_root);
                catch ME
                    warning('Could not create pair folder %s: %s', pair_root, ME.message);
                    continue;
                end
            end

            fprintf('    ↻ WILL RUN: pair [%s,%s] → %s\n', ...
                    J1_str, J2_str, pair_root);

            % Launch the dist/height/tilt sweep for this curvature pair.
            activate_and_monitor_parallel( ...
                template_dir, pair_root, ...
                J1, J2, Prot_rad_var, ...
                dist_vec, height_vec, max_tilt_angle);
        end
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




function activate_and_monitor_parallel(template_dir, output_dir, ...
    J1_curv, J2_curv, Prot_rad_var, dist_vec, height_vec, max_tilt_angle)
% ACTIVATE_AND_MONITOR_PARALLEL
% -------------------------------------------------------------------------
% Wrapper that:
%   • Ensures a parallel pool exists
%   • Opens a diary log
%   • Calls iterative_evolver_run for ONE curvature pair [J1_curv, J2_curv]
%   • Reports total wall-clock time
%
% INPUTS
%   template_dir   – folder with .fe template(s)
%   output_dir     – root folder for this curvature pair (already pair-specific)
%   J1_curv        – curvature of protein type 1 (scalar)
%   J2_curv        – curvature of protein type 2 (scalar)
%   Prot_rad_var   – protein radius (passed downstream if needed)
%   dist_vec       – vector of distances
%   height_vec     – vector of heights
%   max_tilt_angle – maximal tilt (the tilt sampling is done inside
%                    iterative_evolver_run as [-max .. +max] step (2*max/20)
% -------------------------------------------------------------------------

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
        % Call the main function for THIS curvature pair and grid
        iterative_evolver_run( ...
            template_dir, output_dir, ...
            J1_curv, J2_curv, Prot_rad_var, ...
            dist_vec, height_vec, max_tilt_angle);
    catch ME
        fprintf('Error in execution: %s\n', ME.message);

        % Display the function and line number where the error occurred
        for k = 1:length(ME.stack)
            fprintf('Error in function: %s (Line %d)\n', ...
                    ME.stack(k).name, ME.stack(k).line);
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
    J1_curv, J2_curv, rad_P, dist_vec, height_vec, max_tilt_angle)

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % - number of proteins kept for possible future use
    current_num_proteins = 2; %#ok<NASGU>

    % --- Tilt grids for both proteins -----------------------------------
    nTilt_B     = 19;  % protein 1 (B)
    nTilt_E     = 19;  % protein 2 (E)
    tilt_vec_B  = linspace(-max_tilt_angle, max_tilt_angle, nTilt_B);
    tilt_vec_E  = linspace(-max_tilt_angle, max_tilt_angle, nTilt_E);

    % --- timeout parameters (per row of 19 futures) ----------------------
    timeoutSec = 3600;   % max allowed seconds per future
    poll       = 0.25;   % seconds between checks

    % ======= HEIGHT OUTER, DIST INNER ===================================
    for iH = 1:numel(height_vec)
        Height_val = height_vec(iH);

        for iDist = 1:numel(dist_vec)
            Dist_val = dist_vec(iDist);

            % Distance-level folder (prefixed with results_)
            dist_folder = fullfile(output_dir, sprintf('results_Dist_%0.3f', Dist_val));
            if ~exist(dist_folder, 'dir')
                mkdir(dist_folder);
            end

            % Height-level folder (prefixed with results_)
            height_folder = fullfile(dist_folder, sprintf('results_Height_%0.3f', Height_val));
            if ~exist(height_folder, 'dir')
                mkdir(height_folder);
            end

            fprintf('Running J1=%.3f J2=%.3f | Dist=%.3f | Height=%.3f\n', ...
                    J1_curv, J2_curv, Dist_val, Height_val);

            % ---------- LOOP OVER tilt_B (ONE ROW AT A TIME) ----------
            pool = gcp();   % ensure pool

            for iB = 1:nTilt_B
                tilt_B_val = tilt_vec_B(iB);

                % File name for this row (single tiltB, all 19 tiltE)
                rowFileName = sprintf( ...
                    'results_Dist_%0.3f_Height_%0.3f_tiltB_%0.3f_curv_[%.3f,%.3f].mat', ...
                    Dist_val, Height_val, tilt_B_val, J1_curv, J2_curv);
                rowFilePath = fullfile(height_folder, rowFileName);

                % If this row MAT already exists → skip (resume behaviour)
                if exist(rowFilePath, 'file')
                    fprintf('  ✔ SKIP Dist=%.3f Height=%.3f tiltB=%.3f (MAT exists)\n', ...
                            Dist_val, Height_val, tilt_B_val);
                    continue;
                end

                fprintf('  ↻ RUN row: Dist=%.3f Height=%.3f tiltB=%.3f (19 tilt_E)\n', ...
                        Dist_val, Height_val, tilt_B_val);

                % Preallocate row results: 1 × nTilt_E
                rowStruct(1, nTilt_E) = buildNaNStruct( ...
                    rad_P, J1_curv, J2_curv, ...
                    Height_val, Dist_val, ...
                    tilt_B_val, tilt_vec_E(1));

                for iE = 1:nTilt_E
                    tilt_E_val = tilt_vec_E(iE);
                    rowStruct(iE) = buildNaNStruct( ...
                        rad_P, J1_curv, J2_curv, ...
                        Height_val, Dist_val, ...
                        tilt_B_val, tilt_E_val);
                end

                % Launch ONLY 19 futures for this tilt_B row
                futs(nTilt_E,1) = parallel.FevalFuture;
                for iE = 1:nTilt_E
                    tilt_E_val = tilt_vec_E(iE);

                    fprintf(['    [DEBUG] Dist %.3f | Height %.3f | TiltB %.4f | TiltE %.4f | ' ...
                             'J1 %.3f | J2 %.3f\n'], ...
                            Dist_val, Height_val, tilt_B_val, tilt_E_val, J1_curv, J2_curv);

                    futs(iE) = parfeval(pool, @run_single_evolver_iteration, 1, ...
                                        Height_val, ...   % Height argument
                                        Dist_val,  ...    % Distance
                                        J1_curv,   ...    % P1_curv
                                        J2_curv,   ...    % P2_3_curv
                                        tilt_B_val, ...   % tilt of Jcap1
                                        tilt_E_val, ...   % tilt of Jcap2
                                        template_dir, ...
                                        height_folder);
                end

                % ---- watchdog + collection for the 19 results ----
                hasStarted = false(1, nTilt_E);
                isClosed   = false(1, nTilt_E);

                while ~all(isClosed)
                    % A) try to fetch a finished future (times out after poll seconds)
                    idxDone    = [];
                    resultData = [];
                    try
                        [idxDone, resultData] = fetchNext(futs, poll);
                    catch
                        % On timeout / error of fetchNext, see if any finished+errored
                        idxDone = find(strcmp({futs.State}, "finished") & ...
                                       ~cellfun(@isempty, {futs.Error}) & ~isClosed, 1);
                    end

                    % B) mark newly-running jobs (start timer)
                    for kk = find(~hasStarted)
                        if futs(kk).State == "running"
                            hasStarted(kk) = true;
                        end
                    end

                    % C) process finished / errored / cancelled jobs
                    if ~isempty(idxDone) && ~isClosed(idxDone)
                        runSec = seconds(futs(idxDone).RunningDuration);

                        if isempty(futs(idxDone).Error) && runSec < timeoutSec
                            % normal successful completion
                            rowStruct(idxDone) = resultData;
                            fprintf(['    ✔ Dist %.3f | Height %.3f | TiltB %.4f | TiltE %.4f ' ...
                                     'finished in %.2fs\n'], ...
                                    Dist_val, Height_val, ...
                                    tilt_B_val, tilt_vec_E(idxDone), runSec);
                        else
                            % timeout or error
                            cancel(futs(idxDone));
                            if runSec >= timeoutSec
                                msg = sprintf('timeout after %.2fs', runSec);
                                fprintf(2, ['    ⏱ Dist %.3f | Height %.3f | TiltB %.4f | TiltE %.4f ' ...
                                            'exceeded %gs (%.2fs) → NaNs\n'], ...
                                        Dist_val, Height_val, ...
                                        tilt_B_val, tilt_vec_E(idxDone), ...
                                        timeoutSec, runSec);
                            else
                                e   = futs(idxDone).Error;
                                msg = sprintf('error: %s', e.message);
                                fprintf(2, ['    ⚠ Dist %.3f | Height %.3f | TiltB %.4f | TiltE %.4f ' ...
                                            'ERROR: %s → NaNs\n'], ...
                                        Dist_val, Height_val, ...
                                        tilt_B_val, tilt_vec_E(idxDone), e.message);
                            end
                            rowStruct(idxDone).error_message = msg;
                        end

                        isClosed(idxDone) = true;
                        drawnow;
                    end

                    % D) cancel any running job that now exceeds timeout
                    for kk = find(~isClosed & hasStarted)
                        runSec = seconds(futs(kk).RunningDuration);
                        if runSec >= timeoutSec
                            cancel(futs(kk));
                            rowStruct(kk).error_message = "timeout";
                            fprintf(2, ['    ⏱ Dist %.3f | Height %.3f | TiltB %.4f | TiltE %.4f ' ...
                                        'exceeded %gs (%.2fs) → NaNs\n'], ...
                                    Dist_val, Height_val, ...
                                    tilt_B_val, tilt_vec_E(kk), ...
                                    timeoutSec, runSec);
                            isClosed(kk) = true;
                        end
                    end

                    % E) brief pause if some slots still pending
                    if any(~hasStarted & ~isClosed)
                        pause(poll);
                    end
                end % while ~all(isClosed)

                % Convert this 1×19 row into a table and SAVE it now
                subTbl = struct2table(rowStruct, 'AsArray', true);
                save(rowFilePath, 'subTbl', '-v7');

                fprintf('   💾 Saved MAT for Dist=%.3f Height=%.3f tiltB=%.3f\n', ...
                        Dist_val, Height_val, tilt_B_val);
            end % for iB
        end % for iDist
    end % for iH

    disp('All iterations complete - MAT tables written.');
end






function S = buildNaNStruct(radP,J1,J2,Height_val,distVal,tilt_B_val,tilt_E_val)
% Same fields as the real result, but:
%   • geometry / input columns are still filled
%   • everything the worker should have produced is NaN

    if nargin < 6
        tilt_B_val = NaN;
    end
    if nargin < 7
        tilt_E_val = NaN;
    end

    S = struct( ...
      'radius_P2',radP, 'radius_P5',radP, ...
      'Curv_P1_fam', NaN, 'Curv_P2_fam', NaN, ...
      'Curv_P2',J1, 'Curv_P5',J2, ...   % J1,J2 are the family curvatures
      'Height_diff',Height_val, ...
      'distance',distVal, ...
      'tilt_angle_B', tilt_B_val, ...
      'tilt_angle_E', tilt_E_val, ...
      'membrane_tension',NaN, ...
      'total_energy',NaN, ...
      'total_Bending_energy',NaN, ...
      'total_Bending_energy_calculated',NaN, ...
      'plane_bending_energy',NaN, ...
      'cap_bending_energy',NaN, ...
      'total_area',NaN, ...
      'protein_caps_area',NaN, ...
      'membrane_area',NaN, ...
      'protein_base_area_B',NaN, ...
      'protein_base_area_E',NaN, ...
      'total_protein_base_area',NaN, ...
      'analytic_total_area', NaN, ...
      'PHI',NaN, ...
      'error_message','' ...
    );
end



function result = run_single_evolver_iteration(Height_diff_val, Distance, P1_curv, P2_3_curv, tilt_B_val, tilt_E_val, template_dir, output_dir)
    % Initialize result
% --- initialize result fields to mirror OUTPUTLINE names exactly ---
result.radius_P2 = NaN;
result.radius_P5 = NaN;

result.Curv_P1_fam = P1_curv;
result.Curv_P2_fam = P2_3_curv;

result.Curv_P2 = P1_curv;      % (keep your original family logic if intentional)
result.Curv_P5 = P2_3_curv;

result.Height_diff  = Height_diff_val;
result.distance     = Distance;
result.tilt_angle_B = tilt_B_val;   % protein 1 tilt
result.tilt_angle_E = tilt_E_val;   % protein 2 tilt
result.membrane_tension = NaN;

result.total_energy = NaN;
result.total_Bending_energy = NaN;
result.total_Bending_energy_calculated = NaN;
result.plane_bending_energy = NaN;   % NOTE: was "plain_Bending_energy"
result.cap_bending_energy = NaN;

result.total_area = NaN;
result.protein_caps_area = NaN;
result.membrane_area = NaN;

result.protein_base_area_B = NaN;
result.protein_base_area_E = NaN;
result.total_protein_base_area = NaN;

result.analytic_total_area = NaN;

result.PHI = NaN;

result.error_message = '';  % keep for diagnostics

    

    % Full path to the template file in the template directory
    fname_template = fullfile(template_dir, '2_protein_config_square_flat_larger_tilt_template.fe');
    % fname_run = fullfile(output_dir, sprintf('run_curvature_configuration_%d_proteins_sphere_curv_%s_Proteins_curv_[%.3f, %.3f]_distance_%d.fe',current_num_proteins,Curvature, P1_curv, P2_3_curv, Distance));  % Save the run file in the specified directory
    fname_run = fullfile(output_dir, ...
        sprintf('run_dist_%0.3f_height_%0.3f_P1_%0.3f_P23_%0.3f_tiltB_%0.3f_tiltE_%0.3f.fe', ...
        Distance, Height_diff_val, P1_curv, P2_3_curv, tilt_B_val, tilt_E_val));


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

    height_B_VAR = Height_diff_val/2;
    height_E_VAR = -Height_diff_val/2;

    % Replace placeholders with the current variable values
    X = char(X.');
    Y = strrep(X, 'DIST_VAR', num2str(Distance));  % Replace distance placeholder
    Y = strrep(Y, 'TILT_B_VAR', num2str(tilt_B_val));   % Replace P_2_3_CURV_VAR with current Curvature
    Y = strrep(Y, 'TILT_E_VAR', num2str(tilt_E_val));   % Replace P_2_3_CURV_VAR with current Curvature
    Y = strrep(Y, 'HEIGHT_B_VAR', num2str(height_B_VAR));   % Replace P_2_3_CURV_VAR with current Curvature
    Y = strrep(Y, 'HEIGHT_E_VAR', num2str(height_E_VAR));   % Replace P_2_3_CURV_VAR with current Curvature
    Y = strrep(Y, 'Curv_P1_fam_VAR', num2str(P1_curv));   % curv P1 for the family
    Y = strrep(Y, 'Curv_P2_fam_VAR', num2str(P2_3_curv));   % curv P2 for the family
    red = '4';
    blue = '1';
        
    % proteins curvature
    Y = strrep(Y, 'P_B_CURVA_VAR', num2str(P1_curv));
    Y = strrep(Y, 'P_E_CURVA_VAR', num2str(P2_3_curv));
    
    % proteins colors
    Y = strrep(Y, 'P_B_COLOR_VAR', red);
    Y = strrep(Y, 'P_E_COLOR_VAR', blue);
        




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
            result.radius_P2 = NaN;
            result.radius_P5 = NaN;

            result.Curv_P1_fam = P1_curv;
            result.Curv_P2_fam = P2_3_curv;
            
            result.Curv_P2   = P1_curv;
            result.Curv_P5   = P2_3_curv;
            
            result.Height_diff  = Height_diff_val;
            result.tilt_angle_B = tilt_B_val;
            result.tilt_angle_E = tilt_E_val;
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
            
            result.protein_base_area_B      = NaN;
            result.protein_base_area_E      = NaN;
            result.total_protein_base_area  = NaN;

            
            result.analytic_total_area = NaN;
            result.PHI = NaN;
            
            result.error_message = 'model converged wrong (no OUTPUTLINE)';

            else

            for n = 1:numel(linenum)
                tmp = TextAsCells{linenum(n)};
                tmp = strrep(tmp, 'OUTPUTLINE ', '');
                eval(tmp);  % Update the corresponding variables
            end
    
         % Now assign the values extracted from Surface Evolver Captured from OUTPUTLINE in Evolver

            result.radius_P2 = radius_P2;
            result.radius_P5 = radius_P5;

            result.Curv_P1_fam = Curv_P1_fam;
            result.Curv_P2_fam = Curv_P2_fam;
            
            result.Curv_P2 = Curv_P2;
            result.Curv_P5 = Curv_P5;
            
            result.Height_diff = height_diff;
            result.distance = distance;
            result.tilt_angle_B = tilt_angle_B;
            result.tilt_angle_E = tilt_angle_E;
            result.membrane_tension = membrane_tension;
            
            result.total_energy = total_energy;
            result.total_Bending_energy = total_Bending_energy;
            result.total_Bending_energy_calculated = total_Bending_energy_calculated;
            result.plane_bending_energy = plane_bending_energy;   % NOTE: renamed from plain_Bending_energy
            result.cap_bending_energy   = cap_bending_energy;
            
            result.total_area = total_area;
            result.protein_caps_area = protein_caps_area;
            result.membrane_area = membrane_area;

            result.protein_base_area_B = protein_base_area_B;
            result.protein_base_area_E = protein_base_area_E;
            result.total_protein_base_area = total_protein_base_area;
            
            result.analytic_total_area = analytic_total_area;
            result.PHI = PHI;

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

