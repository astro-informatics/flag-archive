function f = flag_synthesis(flmn, L, N, varargin)

% flag_synthesis - Compute Fourier-Laguerre Synthesis
%
% Default usage :
%
%   f = flag_synthesis(flmn, L, N, <options>)
%
% where L and N are the harmonic band-limits, 
% flmn is a complex array of size N x L^2
% The output f is a real or complex array of size N x L*(2*L-1)
% Sampling scheme for theta/phi : McEwen & Wiaux (2011)
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('flmn', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addParamValue('Nodes', 0.0, @isnumeric);
p.addParamValue('Reality', false, @islogical);
p.parse(flmn, L, N, varargin{:});
args = p.Results;

% Compute inverse transform.
f_vec = flag_synthesis_mex(flmn, L, N, args.Nodes, args.Reality);

f = zeros(N, L, (2*L-1));
for n = 1:N
    temp = f_vec(n,:);
    f(n,:,:) = flag_mw_vec2arr( temp );
end
    
end