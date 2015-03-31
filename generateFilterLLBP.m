function filtLLBP = generateFilterLLBP(nElems)
% generate filter suitable for line LBP
% For illustartion see: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3251987/figure/f6-sensors-11-11357/ 
iCenter = round(nElems/2);
nFiltElems = iCenter -1;
filtLLBP = zeros(nElems, 1, nFiltElems);
for iElem = 1:nFiltElems
    %iEnabledElems = [iCenter-iElem, iCenter+iElem];
    iEnabledElems = iCenter-iElem;
    filtLLBP(iEnabledElems, 1, iElem) = 1;
    filtLLBP(iCenter, 1, iElem) = -1;
end