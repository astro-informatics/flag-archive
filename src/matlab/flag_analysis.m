function flmn = flag_analysis(f, L, N, varargin)

% flag_analysis - Compute Fourier-Laguerre Analysis
%
% Default usage :
%
%   flmn = flag_analysis(f, L, N, <options>)
%
% where L and N are the harmonic band-limits, 
% f is a real or complex array of size N x L*(2*L-1)
% The output flmn is a complex array of size N x L^2
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
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addParamValue('Reality', false, @islogical);
p.parse(f, L, N, varargin{:});
args = p.Results;

f_vec = zeros(N, L*(2*L-1));
for n = 1:N
    temp(:,:) = f(n,:,:);
    f_vec(n,:) = flag_mw_arr2vec( temp );
end

flmn = flag_analysis_mex(f_vec, L, N, args.Reality);
  
end