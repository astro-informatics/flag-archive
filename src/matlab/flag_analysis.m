function flmn = flag_analysis(f, varargin)

% flag_analysis - Compute Fourier-Laguerre Analysis
%
% Default usage :
%
%   flmn = flag_analysis(f, <options>)
%
% where L and N are the harmonic band-limits, 
% f is a real or complex array of size N x L*(2*L-1)
% The output flmn is a complex array of size N x L^2
% Sampling scheme for theta/phi : McEwen & Wiaux (2011)
%
% Options :
%  'Reality'         = { false        [do not assume f real (default)],
%                        true         [assume f real (improves performance)] }
%  'L'               = { Harmonic band-limit; L > 1 (default=guessed) }
%  'N'               = { Radial band-limit; N > 1 (default=guessed) }
%  'R'               = { Radial boundary; R > 0 (default=1.0) }
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

sz = size(f);
Nguessed = sz(1);
Lguessed = sz(2);

p = inputParser;
p.addRequired('f', @isnumeric);          
p.addParamValue('L', Lguessed, @isnumeric);          
p.addParamValue('N', Nguessed, @isnumeric);   
p.addParamValue('R', 1.0, @isnumeric);   
p.addParamValue('Reality', false, @islogical);
p.parse(f, varargin{:});
args = p.Results;

f_vec = zeros(args.N, args.L*(2*args.L-1));
for n = 1:args.N
    temp(:,:) = f(n,:,:);
    f_vec(n,:) = flag_mw_arr2vec( temp );
end

flmn = flag_analysis_mex(f_vec, args.L, args.L, args.R, args.Reality);
  
end