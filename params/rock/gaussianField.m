function p = gaussianField(N, vals, sz, std)
%Compute realization of lognormal, isotropic permeability field.
%
% SYNOPSIS:
%   p  = gaussianField(N)
%   p  = gaussianField(N, vals)
%   p  = gaussianField(N, vals, size)
%   p  = gaussianField(N, vals, size, std)
%
% DESCRIPTION:
%   This function creates an approximate Gaussian random field with
%   dimensions given by N. The field is generated by convolving a normal
%   distributed random field with a Gaussian filter.
%
% PARAMETERS:
%   N    - Three-element vector, `[nx, ny, nz]`, specifying the number of
%          discrete values in the `x`, `y`, and `z` coordinate directions
%          respectively.
%
%   VALS - Target interval for the values of the random field
%
%   SIZE - Sets the size of the convolution kernel. Default value is
%          `[3,3,3]`. If `SIZE` is a scalar,  the size is interpreted as 
%          `[SIZE SIZE SIZE].` The convolution kernel must be specified so
%          that numel(SIZE)>=numel(N).
%
%   STD  - Standard deviation used in the Gaussian filter (default: 0.65)
%
% RETURNS:
%   p - The scalar `nx` by `ny` by `nz` Gaussian random field

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


if nargin==1
   sz = 3; std = 0.65; vals = [-1 1];
elseif nargin==2
   sz = 3; std = 0.65;
elseif nargin==3
   std = 0.65;
elseif nargin>4 || nargin==0
   error('gaussianField:wrongNumberOfInputs','Wrong number of input arguments.');
end
if length(sz)==1
   sz = [sz sz sz];
elseif numel(sz) < numel(N)
    error('gaussianField:invalidSizeInput','SIZE must have the same number of components as N');
elseif numel(sz)>3 || numel(sz)==0
   error('gaussianField:invalidSizeInput','SIZE must be a vector with 1, 2, or 3 elements.');
end
padSz = (sz-1)/2;

% random start vector
p = randn(N);

% construct Gaussian filter
if numel(N)==1
   x = -padSz(1):padSz(1);
   D = exp(-x.*x/(2*std*std));
elseif numel(N)==2
   [x,y] = meshgrid(-padSz(2):padSz(2),-padSz(1):padSz(1));
   D = exp(-(x.*x + y.*y)/(2*std*std));
elseif numel(N)==3
   [x,y,z] = meshgrid(-padSz(2):padSz(2),-padSz(1):padSz(1),-padSz(3):padSz(3));
   D = exp(-(x.*x + y.*y + z.*z)/(2*std*std));
else
   error(id('InvalidDimInput'),'N must be a vector with 1, 2, or 3 elements.');
end
D = D / sum(D(:));

% Pad geometry before convolution
idx = cell(numel(N),1);
for k=1:numel(N)
   M = size(p,k); I = ones(1,padSz(k)); idx{k} = [I 1:M M*I];
end

% Convolve filter with data
p  = convn(p(idx{:}), D, 'valid');

% Scale to correct interval
if numel(vals)~=2
   error(id('InvalidValsInput'),'VALS must be a 2 element vector');
end
mp = min(p(:)); Mp = max(p(:));
p  = vals(1)  + (vals(2)-vals(1))*(p - mp)/(Mp-mp);
