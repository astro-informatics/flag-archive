function fn = slag_analysis(f, varargin)

% slag_analysis - Compute 1D spherical Laguerre Analysis
%
% Default usage :
%
%   fn = slag_analysis(f, P, tau)
%
% where P is the harmonic band-limit, 
% tau is the radial scale factor,
% f is a vector of size P,
% the output fn is a vector of P
%
% Options :
%  'P'    = { Band-limit; P > 1 (default=guessed) }
%  'tau'    = { double (default=1.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Pguessed = max([sz(1) sz(2)]);

p = inputParser;
p.addRequired('f', @isnumeric);          
p.addParamValue('P', Pguessed, @isnumeric);   
p.addParamValue('tau', 1.0, @isnumeric);  
p.parse(f, varargin{:});
args = p.Results;

fn = slag_analysis_mex(f, args.P, args.tau);
    
end