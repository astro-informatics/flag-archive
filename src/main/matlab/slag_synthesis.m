function [f, nodes] = slag_synthesis(fn, varargin)

% slag_synthesis - Compute 1D spherical Laguerre Synthesis
%
% Default usage :
%
%   f = slag_synthesis(fn, P, <options>)
%
% where P is the harmonic band-limit, 
% fn is a vector of size P,
% the output f is a vector of size P,
% You must specify either the noded of the sampling
% or the radial scale factor tau (default: 1.0)
%
% Options :
%  'P'    = { Band-limit; P > 1 (default=guessed) }
%  'tau'    = { double (default=1.0) }
%  'nodes'    = { double (default=0.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(fn);
Pguessed = max([sz(1) sz(2)]);

p = inputParser;
p.addRequired('fn', @isnumeric);          
p.addParamValue('P', Pguessed, @isnumeric);   
p.addParamValue('Nodes', 0.0, @isnumeric);
p.addParamValue('tau', 1.0, @isnumeric);
p.parse(fn, varargin{:});
args = p.Results;

[f, nodes] = slag_synthesis_mex(fn, args.P, args.tau, args.Nodes);
  
end