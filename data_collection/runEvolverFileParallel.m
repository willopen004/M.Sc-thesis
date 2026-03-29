function [status, out_string] = runEvolverFileParallel(fe_name, varargin)
    %                           = runEvolverFile(fe_name, run_command, fe_path, evolver_path, interactive, quiet_load)
    % Inputs:
    % fe_name = evolver script file to load
    % run_command = evolver script command string to execute after load
    % fe_path = path of evolver file
    % evolver_path = path of evolver program
    % interactive = run surface evolver in interactive mode
    % quiet_load = suppress display of evolver file while loading
    %
    % Outputs:
    % status = status returned from the OS. status=0 indicates evolver executable was run successfully
    % out_string = output string from the evolver (whatever appears on screen when running evolver as stand-alone)

    % Set default values
    run_command = 'runme();';
    fe_path = '';
    evolver_path = pwd;
    quiet_load = true;
    interactive = false;

    % Handle varargin safely
    num_args = length(varargin);
    if num_args >= 1 && ~isempty(varargin{1})
        run_command = varargin{1};  % Assign run_command from varargin
    end
    if num_args >= 2 && ~isempty(varargin{2})
        fe_path = varargin{2};  % Assign fe_path from varargin
    end
    if num_args >= 3 && ~isempty(varargin{3})
        evolver_path = varargin{3};  % Assign evolver_path from varargin
    end
    if num_args >= 4 && ~isempty(varargin{4})
        interactive = varargin{4};  % Assign interactive mode from varargin
    end
    if num_args >= 5 && ~isempty(varargin{5})
        quiet_load = varargin{5};  % Assign quiet_load flag from varargin
    end

    % Ensure correct path format for Windows
    if ispc
        fe_path = strrep(fe_path, '/', '\');
    end

    % Setup the command execution string for Surface Evolver
    if interactive
        if isunix
            command_string = [' -x -r "' run_command ';exit 0;" '];
        else
            command_string = [' -r "' run_command ';" '];
        end
    else
        command_string = [' -x -r "' run_command ';exit 0;" '];
    end

    if quiet_load
        command_string = [' -Q ' command_string];
    end

    % Setup the system command string
    evolver_exe = '\evolver.exe';  % Adjust for your evolver executable path
    system_string = ['"' evolver_path evolver_exe ' ' command_string ' ' fe_path fe_name '"'];

    % Create a .cmd file instead of directly running the command
    % **Generate a guaranteed unique .cmd filename**
    cmd_filename = tempname(pwd);  % Generates a unique name in the current directory
    cmd_filename = [cmd_filename, '.cmd'];  % Append .cmd extension

    % >>> ADD – begin  ---------------------------------------------------------
    cleanupObj = onCleanup(@() safeDelete(cmd_filename));
    % >>> ADD – end  

    fid = fopen(cmd_filename, 'wt');  % Open the file for writing

    if fid == -1
        error('Failed to create unique Surface Evolver command file.');
    end

    % Write the command into the .cmd file
    fprintf(fid, '@echo off\n');
    fprintf(fid, '"%s%s" %s "%s%s"\n', evolver_path, evolver_exe, command_string, fe_path, fe_name);
    fclose(fid);

    % Run the .cmd file and capture the output
    [status, out_string] = system(cmd_filename);

    % **Delete the temporary .cmd file after execution**
    delete(cmd_filename);

    if status == 0
        disp('Surface Evolver ran successfully.');
    else
        disp('There was an error running Surface Evolver.');
    end
end

function safeDelete(fname)
    if exist(fname,'file'), try delete(fname); catch, end, end
end

