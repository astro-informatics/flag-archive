function fn = slag_analysis(f, N, R)

% slag_analysis - Compute 1D spherical Laguerre Analysis
%
% Default usage :
%
%   fn = slag_analysis(f, N, R)
%
% where N is the harmonic band-limit, 
% R is the radial limit,
% f is a vector of size N,
% the output fn is a vector of N
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('f', @isnumeric);          
p.addRequired('N', @isnumeric);
p.addRequired('R', @isnumeric);   
p.parse(f, N, R);
args = p.Results;

fn = slag_analysis_mex(f, N, R);
  