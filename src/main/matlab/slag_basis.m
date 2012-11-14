function f = slag_basis(N, nodes, tau)

% slag_basis_mex - Compute spherical Laguerre basis function
% of order N on a grid of values nodes (with rescaling factor tau)
%
% Default usage :
%
%   f = slag_basis(N, nodes, tau)
%
% where N is the order of the spherical Laguerre mode,
% nodes is the radial grig on which to compute the basis function,
% tau is the rescaling factor.
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

f = slag_basis_mex(N, nodes, tau);
  
end