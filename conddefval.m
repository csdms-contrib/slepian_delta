% conddefval - Conditionally define a value
%
% This function takes an input value and a default value. If the input value
% is empty, it returns the default value. Otherwise, it returns the input value.
%
% Syntax:
%   value = conddefval(inputValue, defaultValue)
%
% Input Arguments:
%   - inputValue: The input value to be checked.
%   - defaultValue: The default value to be returned if the input value is empty.
%
% Output Argument:
%   - value: The resulting value based on the condition.
%
% Example:
%   value = conddefval([], 10);
%   % Returns 10
%
%   value = conddefval(5, 10);
%   % Returns 5
%
% Author: 'Will' En-Chi Lee
% Last Modified: 2024-05-30
%
function value = conddefval(inputValue, defaultValue)

    if isempty(inputValue)
        value = defaultValue;
    else
        value = inputValue;
    end

end
