function f = slag_basis(N, nodes, tau)

% slag_basis_mex - Compute the spherical Laguerre basis functions
% up to order N on a grid of radii nodes (with rescaling factor tau)
%
% Default usage :
%
%   f = slag_basis(N, nodes, tau)
%   plot(f)
%
% where N is the maximum order/mode of the spherical Laguerre functions,
% nodes is the radial grig on which to compute the basis function,
% tau is the rescaling factor.
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

f = slag_basis_mex(N, nodes, tau);
  
end