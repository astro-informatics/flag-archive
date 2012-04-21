function [nodes, weights] = slag_gausslaguerre_quadrature(N, alpha)

% slag_syslag_gausslaguerre_quadraturenthesis - 
% Compute Gauss-Laguerre quadrature nodes and weights
%
% Default usage :
%
%   [nodes, weights] = slag_gausslaguerre_quadrature(N, alpha)
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

[nodes, weights] = slag_gausslaguerre_quadrature_mex(N, alpha);
  
end