

function variable_list = add_var(variable_list, new_variable)
% variable_list = add_var(variable_list, new_variable)
%
% PURPOSE: 
% add new variable name to the list of current save variables
% so that they do not repeat
%
% INPUTS
%       'variable_list' : cell array of variables to save and not erase
%       'new_variable' : character variable with name of variable to save
%               e.g. 'data'
%
% OUTPUTS
%       'variable_list' : updated variable list with the new variable (if
%               not already present)
%
% EXAMPLE USE
%       initial_var = add_var(initial_var, 'data')
%       clearvars('-except',initial_var{:})
%
% ES DICKINSON 2026


%%

if ~any(strcmp(new_variable, variable_list))
    variable_list{end+1} = new_variable;
end
