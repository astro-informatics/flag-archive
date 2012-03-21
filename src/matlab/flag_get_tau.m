function tau = flag_get_tau(N, R)

% flag_get_tau - Compute scaling factor of SLAG transform
%
% Default usage :
%
%   tau = flag_get_tau(N, R)
%
% where N is the radial harmonic band-limit, 
% R is the radial limit,
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('N', @isnumeric);   
p.addRequired('R', @isnumeric); 
p.parse(N, R);
args = p.Results;

tau = flag_get_tau_mex(N, R);
       
end