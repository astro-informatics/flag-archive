function [rs, thetas, phis] = flag_sampling(L, N, R)

% flag_sampling - Compute Fourier-Laguerre sampling scheme
%
% Default usage :
%
%   [rs, thetas, phis] = flag_sampling(L, N, R)
%
% where L and N are the harmonic band-limits, 
% R is the radial limit,
% Output :
%   rs is a real vector of size N
%   thetas is a real vector of size L
%   phis is a real vector of size 2*L-1
% Sampling scheme for theta/phi : McEwen & Wiaux (2011)
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('L', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addRequired('R', @isnumeric); 
p.parse(L, N, R);
args = p.Results;

[rs, thetas, phis] = flag_sampling_mex(L, N, R);
     