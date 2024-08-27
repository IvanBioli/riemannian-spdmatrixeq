% Description: Adds subfolders to MATLAB's path
run manopt/importmanopt.m
addpath(genpath("baseline_solvers"))
addpath(genpath("examples"))
addpath(genpath("manifold_optimization"))
addpath(genpath("paper_experiments"))
addpath(genpath("utils"))
fprintf('The folders needed to run the code were added to Matlab''s path.\n');
response = input('Save path for future Matlab sessions? [Y/N] ', 's');
if strcmpi(response, 'Y')
    failed = savepath();
    if ~failed
        fprintf('Path saved: no need to call init.m next time.\n');
    else
        fprintf(['Something went wrong.. Perhaps missing permission ' ...
            'to write on pathdef.m?\nPath not saved: ' ...
            'please re-call init.m next time.\n']);
    end
else
    fprintf('Path not saved: please re-call init.m next time.\n');
end