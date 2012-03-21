function [f, nodes] = slag_synthesis(fn, N, varargin)

% slag_synthesis - Compute 1D spherical Laguerre Synthesis
%
% Default usage :
%
%   f = slag_synthesis(fn, N, <options>)
%
% where N is the harmonic band-limit, 
% fn is a vector of size N,
% the output f is a vector of size N,
% You must specify either the noded of the sampling
% or the radial limit R (default: 1.0)
%
% Option :
%  'R'    = { double (default=1.0) }
%  'nodes'    = { double (default=0.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('fn', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addParamValue('Nodes', 0.0, @isnumeric);
p.addParamValue('R', 1.0, @isnumeric);
p.parse(fn, N, varargin{:});
args = p.Results;

[f, nodes] = slag_synthesis_mex(fn, N, args.R, args.Nodes);
  
end