function nodes = slag_sampling(N, R)

% slag_sampling - Compute 1D spherical Laguerre Sampling scheme
%
% Default usage :
%
%   nodes = slag_sampling(N, R)
%
% where N is the harmonic band-limit, 
% R is the radial limit,
% the output nodes is a vector of size N
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('N', @isnumeric);  
p.addRequired('R', @isnumeric);   
p.parse(N, R);
args = p.Results;

nodes = slag_sampling_mex(N, R);
  