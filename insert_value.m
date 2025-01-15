function [adjusted_vector] = insert_value(vector, val)
%% DESCRIPTION: This function will insert a value into a vector at the front and remove from the back
%
%% INPUT:
% vector                    The vector to adjust
%
%% OUTPUT: 
% adjusted_vector           The adjusted vector
%% Code
adjusted_vector = [val];
for i = 1:size(vector, 2)-1
    adjusted_vector = [adjusted_vector vector(:, i)];
end

