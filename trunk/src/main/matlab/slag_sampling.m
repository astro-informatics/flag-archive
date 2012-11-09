function nodes = slag_sampling(P, R)

% slag_sampling - Compute 1D spherical Laguerre Sampling scheme
%
% Default usage :
%
%   nodes = slag_sampling(P, R)
%
% where P is the harmonic band-limit, 
% R is the radial limit,
% the output nodes is a vector of size P
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('P', @isnumeric);  
p.addRequired('R', @isnumeric);   
p.parse(P, R);
args = p.Results;

nodes = slag_sampling_mex(P, R);
    
end