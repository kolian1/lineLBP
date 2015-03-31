function filtLLBP = generateShiftLLBP(nElems, isVer)
if nargin<2 % default is Horizontal filter
    isVer = false;
end
% generate filter suitable for line LBP
% For illustartion see: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3251987/figure/f6-sensors-11-11357/ 

% Rows- shift in x- horizontal (left-right) 
% Columns- shift in y- vertical (up-down)

filtLLBP = zeros(nElems, 2);
filtLLBP(1:nElems, 1) = 1:nElems;

if isVer % shift dimentsions
    filtLLBP = filtLLBP(:, [2, 1]);
end