function f = flag_synthesis(flmn, varargin)

% flag_synthesis - Compute Fourier-Laguerre Synthesis
%
% Default usage :
%
%   f = flag_synthesis(flmn, <options>)
%
% where L and N are the harmonic band-limits, 
% flmn is a complex array of size N x L^2
% The output f is a real or complex array of size N x L*(2*L-1)
% Sampling scheme for theta/phi : McEwen & Wiaux (2011)
%
% Option :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed) }
%  'N'               = { Radial band-limit; N > 1 (default=guessed) }
%  'R'               = { Radial boundary; R > 0 (default=1.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(flmn);
Nguessed = sz(1);
Lguessed = sqrt(sz(2));

p = inputParser;
p.addRequired('flmn', @isnumeric);          
p.addParamValue('L', Lguessed, @isnumeric);          
p.addParamValue('N', Nguessed, @isnumeric);   
p.addParamValue('Nodes', 0.0, @isnumeric);
p.addParamValue('R', 1.0, @isnumeric);   
p.addParamValue('Reality', false, @islogical);
p.parse(flmn, varargin{:});
args = p.Results;

% Compute inverse transform.
f_vec = flag_synthesis_mex(flmn, args.L, args.N, args.Nodes, args.R, args.Reality);

sz = size(args.Nodes);
P = sz(2);
if P > 1 
    dim = P ;
else
    dim = args.N + 1 ;
end
    
f = zeros(dim, args.L, (2*args.L-1));
for n = 1:dim
    temp = f_vec(n,:);
    f(n,:,:) = flag_mw_vec2arr( temp );
end
    
end