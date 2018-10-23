function [] = dispcell(A);
% function [] = dispcell(A);
%
% displays elements of a "row" or "column" cell array, A.  Cell array
% elements must be a horizonal array of characters or numbers.
%

s = size(A);
if ~isequal(prod(s),max(s))
   error('only "column" cell array are valid')
end

for j = 1:max(s)
   if ischar(A{j})
      disp(A{j})
   else
      disp(num2str(A{j}))
   end
end
