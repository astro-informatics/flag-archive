function f = spherbessel_basis(ell, nodes)

% spherbessel_basis - Compute spherical Bessel basis function
% of order ell on a grid of values nodes.
%
% Default usage :
%
%   f = spherbessel_basis(ell, nodes)
%
% where N is the order of the spherical Bessel mode,
% nodes is the radial grig on which to compute the basis function,
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

f = spherbessel_basis_mex(ell, nodes);
  
end