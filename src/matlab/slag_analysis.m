function fn = slag_analysis(f, varargin)

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
% Options :
%  'N'    = { Band-limit; N > 1 (default=guessed) }
%  'R'    = { double (default=1.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Nguessed = max([sz(1) sz(2)]) - 1;

p = inputParser;
p.addRequired('f', @isnumeric);          
p.addParamValue('N', Nguessed, @isnumeric);   
p.addParamValue('R', 1.0, @isnumeric);  
p.parse(f, varargin{:});
args = p.Results;

fn = slag_analysis_mex(f, args.N, args.R);
    
end