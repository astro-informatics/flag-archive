function tau = flag_get_tau(P, R)

% flag_get_tau - Compute scaling factor of SLAG transform
%
% Default usage :
%
%   tau = flag_get_tau(P, R)
%
% where P is the radial harmonic band-limit, 
% R is the radial limit,
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

p = inputParser;
p.addRequired('P', @isnumeric);   
p.addRequired('R', @isnumeric); 
p.parse(P, R);
args = p.Results;

tau = flag_get_tau_mex(P, R);
       
end