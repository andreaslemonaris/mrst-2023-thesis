function [info, present] = dataset_spe1()
% Info function for SPE1 dataset. Use getDatasetInfo or getAvailableDatasets for practical purposes.

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    [info, present] = datasetInfoStruct(...
        'name', 'SPE1', ...
        'website', '', ...
        'fileurl', 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/SPE1.zip', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'cells',   300, ...
        'examples', {'ad-blackoil:runAndCompareSPE1', ...
                     'ad-blackoil:demonstrateExtraOutputSPE1', ...
                     'ad-blackoil:timestepComparisonSPE1'}, ...
        'description', 'This is the first SPE benchmark. It is a live oil/dry gas black-oil model with nearly immobile water. There is a single bhp-controlled producer, and a single injector injecting gas into the reservoir initially filled with undersaturated oil.',...
        'filesize',    0.005, ...
        'modelType', 'Three-phase, black-oil with gas in oil phase. Cartesian grid' ...
         );
end
