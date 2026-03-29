
Curvature_model = 'Intrinsic Curvature';   %  'Intrinsic Curvature' | 'Spontaneous Curvature'
T = concatenateMatFilesRows(Curvature_model);   % T is the combined table

function resultTable = concatenateMatFilesRows(Curvature_model)
%CONCATENATEMATFILESROWS  Recursively load every results*.mat under the
% given root folder, harmonise the columns and concatenate all rows into
% one big table.  Saves a copy as .mat and .csv next to the data.
%
%  • Any non-numeric column is converted with STR2DOUBLE.
%  • Columns that are *all* NaN are dropped.
%  • Files with fewer columns than the first “reference” file are skipped.
%
%  Yiftach Navot · Aug-2025
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 1)  Resolve base directory
% -------------------------------------------------------------------------
if strcmpi(Curvature_model,'Intrinsic Curvature')
    baseDirectoryPath = sprintf('E:/Yiftach/OneDrive - Technion/kaiser/my_output_directory_parallel/half_square/V1/');
elseif strcmpi(Curvature_model,'Spontaneous Curvature')
baseDirectoryPath = sprintf('D:/yiftach_OneDrive/OneDrive - Technion/kaiser/my_output_directory_parallel/Intrinsic_Curvature/same_Jcap_V7_fixed_model/');
else
    error('Unknown Curvature_model "%s"',Curvature_model);
end
baseDirectoryPath = char(baseDirectoryPath);   % ensure char for DIR

% -------------------------------------------------------------------------
% 2)  Locate *.mat files (depth-first, any sub-folder)
% -------------------------------------------------------------------------
matFiles = dir(fullfile(baseDirectoryPath,'**','results*.mat'));
if isempty(matFiles)
    error('No results*.mat files found under "%s"',baseDirectoryPath);
end

% -------------------------------------------------------------------------
% 3)  Loop through files
% -------------------------------------------------------------------------
resultTable   = table();  % grows as we append
commonHeaders = {};

for f = 1:numel(matFiles)
    filePath = fullfile(matFiles(f).folder,matFiles(f).name);

    % ---- 3.a  load ONE table variable from that .mat --------------------
    S = load(filePath);                         % struct of vars
    tblLocal = [];
    for vn = fieldnames(S).'
        if istable(S.(vn{1}))
            tblLocal = S.(vn{1});
            break
        end
    end
    if isempty(tblLocal)
        warning('Skipping %s – no table variable inside.',filePath);
        continue
    end

    % ---- 3.b  clean table (helper below) --------------------------------
    if f==1
        [tblLocal, commonHeaders] = cleanTable(tblLocal, [], filePath);  % first file
        resultTable             = tblLocal;
        if ~ismember('error_message', commonHeaders)
            % add an empty string column and register it in the header list
            tblLocal.error_message = strings(height(tblLocal),1);
            commonHeaders{end+1}    = 'error_message';
        end
    else
        tblLocal = cleanTable(tblLocal, commonHeaders, filePath);    % later files
        if isempty(tblLocal)                               % header mis-match
            continue, end
        resultTable = [resultTable ; tblLocal];            %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
% 4)  Drop rows that are entirely NaN (all numeric cols missing)
% -------------------------------------------------------------------------
allNanRow = all(ismissing(resultTable),2);
resultTable(allNanRow,:) = [];

% -------------------------------------------------------------------------
% 5)  Save as MAT and CSV next to the original data
% -------------------------------------------------------------------------
timestamp   = datestr(now,'dd.mm.yyyy_HH.MM');
matName     = sprintf('data_table_%s_%s.mat',   replace(Curvature_model,' ','_'),timestamp);
csvName     = sprintf('data_table_%s_%s.csv',   replace(Curvature_model,' ','_'),timestamp);
matFullPath = fullfile(baseDirectoryPath,matName);
csvFullPath = fullfile(baseDirectoryPath,csvName);

save(matFullPath,'resultTable','-v7');
writetable(resultTable,csvFullPath);

fprintf('\n✅  Combined table written:\n   %s\n',matFullPath);
fprintf('   %s\n',csvFullPath);
end




function [tbl, headers] = cleanTable(tbl, referenceHeaders, filePath)
% ------------------------------------------------------------------
% if this is the FIRST file  →  no second argument was supplied
% ------------------------------------------------------------------
if nargin < 3
    filePath = '(unknown file)';
end

% -------- detect first call -------------------------------------------
isFirst = (nargin < 2) || isempty(referenceHeaders);

if isFirst
    % ---------- FIRST-FILE LOGIC (unchanged) ---------------------------
    % drop Var… columns
    bad = startsWith(tbl.Properties.VariableNames,'Var');
    tbl(:,bad) = [];

    % ensure error_message column exists and is string
    if ~ismember('error_message',tbl.Properties.VariableNames)
        tbl.error_message = strings(height(tbl),1);
    else
        tbl.error_message = string(tbl.error_message);
    end

    % convert non-numeric columns to double
    for k = 1:width(tbl)
        if ~strcmp(tbl.Properties.VariableNames{k},'error_message') ...
                && ~isnumeric(tbl.(k))
            tbl.(k) = str2double(string(tbl.(k)));
        end
    end

    % drop all-NaN columns (except error_message)
    nanCols = all(ismissing(tbl),1) & ...
              ~strcmp(tbl.Properties.VariableNames,'error_message');
    tbl(:,nanCols) = [];

    headers = tbl.Properties.VariableNames;
    return                     % ---- finished first file -------------
end
% ------------------------------------------------------------------
% If we get here we *do* have referenceHeaders: second, third, … file
% ------------------------------------------------------------------
headers = referenceHeaders;          % keep output count consistent
colNames = tbl.Properties.VariableNames;

% same “Var…” / NaN / type-forcing code as above -------------------
bad = startsWith(colNames,'Var');   tbl(:,bad) = [];

for k = 1:width(tbl)
    if strcmp(colNames{k},'error_message')
        tbl.(k) = string(tbl.(k));
    elseif ~isnumeric(tbl.(k))
        tbl.(k) = str2double(string(tbl.(k)));
    end
end
nanCols = all(ismissing(tbl),1) & ~strcmp(colNames,'error_message');
tbl(:,nanCols) = [];

% ----- add missing / drop extras / enforce order -------------------
miss = setdiff(referenceHeaders,colNames);
for m = miss
    if strcmp(m{1},'error_message')
        tbl.(m{1}) = strings(height(tbl),1);
    else
        tbl.(m{1}) = NaN(height(tbl),1);
    end
end
extras = setdiff(tbl.Properties.VariableNames,referenceHeaders);
tbl(:,extras) = [];

% ───── reorder columns to the reference order ─────────────
try
    tbl = tbl(:, referenceHeaders);   % <-- possible failure point
catch
    missingCols = setdiff(referenceHeaders, tbl.Properties.VariableNames);
    extraCols   = setdiff(tbl.Properties.VariableNames, referenceHeaders);
    warning('⚠ Skipping file: %s\n   • Missing columns: %s\n   • Extra columns: %s', ...
            filePath, strjoin(missingCols, ', '), strjoin(extraCols, ', '));
    tbl = table();     % return empty → caller will skip this file
end

end

