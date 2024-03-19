function mult = accumulateCartesianMultipliers(G, mult)
%Compute Face Multipliers From Collection of Cell/Face Multipliers
%
% SYNOPSIS:
%   faceMult = accumulateCartesianMultipliers(G, cellFaceMult)
%
% PARAMETERS:
%   G            - Grid structure.  Must identify Cardinal/Cartesian
%                  direction in column 2 of `G.cells.faces`.  Typically
%                  created by construction function `processGRDECL` or
%                  `tensorGrid`.
%
%   cellFaceMult - Collection of multipliers defined in terms of pairs of
%                  cells and Cardinal directions.  Structure with one or
%                  more of the following fields, each assumed to be numeric
%                  arrays of `G.cells.num` or `PROD(G.cartDims)` elements:
%
%                    - x and x_: Multipliers per cell in cell's positive
%                      and negative X directions respectively.
%
%                    - y and y_: Multipliers per cell in cell's positive
%                      and negative Y directions respectively.
%
%                    - z and z_: Multipliers per cell in cell's positive
%                      and negative Z directions respectively.
%
%                  This multiplier structure often corresponds to the
%                  `rock.multipliers` field generated by `grdecl2Rock` or
%                  similar functions.
%
% RETURNS:
%   faceMult - Numeric array, size `G.faces.num`-by-1, of aggregate
%              multiplier values.  In particular, the `i`-th element is the
%              product of all multiplier values from the input data that
%              apply to grid face `i`.  This value is guaranteed to be one
%              (1) if no input multiplier applies to the `i`-th grid face.
%
%              If there either are no multipliers in the input array at all
%              or if the input grid does not identify Cartesian directions
%              then faceMult is an empty array of size `G.faces.num`-by-0.
%
% EXAMPLE:
%   % Input 'NORNE' dataset and build a grid and a rock structure for it.
%   % This structure contains transmissibility multipliers in the positive
%   % Z direction (MULTZ, downwards).
%   %
%   grdecl = readGRDECL(fullfile(getDatasetPath('NORNE'), 'NORNE.GRDECL'));
%   grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
%   G      = processGRDECL(grdecl, 'SplitDisconnected', false);
%   rock   = grdecl2Rock(grdecl, G.cells.indexMap);
%
%   % Aggregate multiplier values to grid's active faces.
%   mult = accumulateCartesianMultipliers(G, rock.multipliers);
%
%   % Visualise layered multiplier structure by colouring those faces where
%   % the total transmissibility multiplier value is less than 0.5 by the
%   % base 10 logarithm of the total multiplier value.
%   %
%   pick = mult < 0.5;
%   plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.25)
%   plotFaces(G, pick, log10(mult(pick)), 'FaceAlpha', 0.625)
%   view(-300, 13), axis tight, grid on
%
% SEE ALSO:
%   `tensorGrid`, `processGRDECL`, `grdecl2Rock`.

%{
Copyright 2020-2023 SINTEF Digital, Mathematics & Cybernetics.

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

   multarr = normalise_multiplier_array(G, mult);

   if isempty(multarr)
      mult = multarr;
      return
   end

   mult = accumulate(G, multarr);
end

%--------------------------------------------------------------------------

function multarr = normalise_multiplier_array(G, mult)
   if isempty(mult) || (size(G.cells.faces, 2) < 2)
      multarr = ones([G.faces.num, 0]);

      if ~isempty(mult)
         warning('CartesianDirections:Missing', ...
                ['Input grid does identify Cartesian cell directions ', ...
                 '(second column of ''cells.faces'' missing)']);
      end
   else
      mfld     = fieldnames(mult);
      id       = identify_directions();
      totcells = prod(G.cartDims);

      multarr.dir  = zeros([1, numel(mfld)]);
      multarr.mult = zeros([G.cells.num, numel(mfld)]);

      for k = 1 : numel(mfld)
         try
            dir = id(mfld{k});
         catch
            warning('Multiplier:Unknown', ...
                   ['Unknown multiplier direction ', ...
                    'field ''%s''.  Ignored.'], mfld{k});
            continue
         end

         value = mult.(mfld{k});

         if numel(value) == G.cells.num
            multarr.dir(k)     = dir;
            multarr.mult(:, k) = value;

         elseif numel(value) == totcells
            multarr.dir(k)     = dir;
            multarr.mult(:, k) = value(G.cells.indexMap);

         else
            warning('Multiplier:IncompatibleSize', ...
                   ['Multiplier array ''%s'' has incompatible size ', ...
                    '''%d''.  Expected one of ''%d'' or ''%d''.\n', ...
                    '  -> Multiplier array ignored.'], ...
                    mfld{k}, numel(value), G.cells.num, totcells);
         end
      end
   end
end

%--------------------------------------------------------------------------

function mult = accumulate(G, m)
   assert ((min(m.dir(:)) >= 1) && (max(m.dir(:)) <= 2 * G.griddim), ...
           'Incompatible cell-face directions');

   % m.dir  is 1-by-k array of cell-face directions in 1 .. 2*griddim.
   % m.mult is normalised G.cells.num-by-k array of input multipliers.

   % Each column of 'pick' has exactly one `true` element, corresponding to
   % the direction of that column/element in m.dir.
   ndir = numel(m.dir);
   pick = false([2 * G.griddim, ndir]);
   pick(sub2ind(size(pick), m.dir(:), (1 : ndir).')) = true;

   cn = gridCellNo(G);
   f  = repmat(G.cells.faces(:,1), [ndir, 1]);
   m  = reshape(m.mult(cn, :), [], 1);

   facedir  = repmat(G.cells.faces(:,2), [ndir, 1]);
   columnid = rldecode((1 : ndir) .', numel(cn));
   i = pick(sub2ind(size(pick), facedir, columnid));

   mult = accumarray(f(i), m(i), [G.faces.num, 1], @prod, 1);
end

%--------------------------------------------------------------------------

function id = identify_directions()
   s  = struct('x_', 1, 'i_', 1, 'x', 2, 'i', 2, ...
               'y_', 3, 'j_', 3, 'y', 4, 'j', 4, ...
               'z_', 5, 'k_', 5, 'z', 6, 'k', 6);
   id = @(d) s.(lower(regexprep(d, '\W', '_')));
end
