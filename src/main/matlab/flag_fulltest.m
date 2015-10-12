% flag_fulltest - tauun all tests
%
% FLAG package to perform 3D Fourier-Laguerre Analysis
% Copyright (C) 2012  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

clear all;
close all;

% Main parameters
L = 32
P = 32
tau = 1.0    

% Generate random 3D FLAG decomposition
flmn = zeros(P, L^2);
flmn = rand(size(flmn)) + sqrt(-1)*rand(size(flmn));
flmn = 2.*(flmn - (1+sqrt(-1))./2);

% Test exactness for default FLAG transform
f = flag_synthesis(flmn);
flmn_rec = flag_analysis(f);
flag_default_transform_error = max(max(abs(flmn-flmn_rec)))

% Test exactness for complex FLAG transform
f = flag_synthesis(flmn, 'L', L, 'P', P);
flmn_rec = flag_analysis(f, 'L', L, 'P', P);
flag_custom_transform_error = max(max(abs(flmn-flmn_rec)))

nodes = slag_sampling(P, tau);
% Test exactness for complex FLAG transform on the same grid
f = flag_synthesis(flmn, 'L', L, 'P', P, 'Nodes', nodes);
flmn_rec = flag_analysis(f, 'L', L, 'P', P, 'tau', tau);
flag_grid_transform_error = max(max(abs(flmn-flmn_rec)))

nodes = slag_sampling(P, tau);
% Test exactness for complex FLAG transform on the same grid
f = flag_synthesis(flmn, 'tau', tau);
flmn_rec = flag_analysis(f, 'tau', tau);
flag_bound_transform_error = max(max(abs(flmn-flmn_rec)))

% Impose reality on flms.
for en = 1:P
   for el = 0:L-1
      ind = el*el + el + 1;
      flmn(en,ind) = real(flmn(en,ind));
      for m = 1:el
         ind_pm = el*el + el + m + 1;
         ind_nm = el*el + el - m + 1;
         flmn(en,ind_nm) = (-1)^m * conj(flmn(en,ind_pm));
      end  
   end
end
% Test exactness of real FLAG transform
f = flag_synthesis(flmn, 'L', L, 'P', P, 'reality', true);
flmn_rec = flag_analysis(f, 'L', L, 'P', P, 'reality', true);
flag_real_transform_error = max(max(abs(flmn-flmn_rec)))

% Generate random 1D SLAG decomposition
fn = rand(1,P);

% Test exactness of SLAG transform
[f, nodes] = slag_synthesis(fn);
fn_rec = slag_analysis(f);
slag_default_transform_error = max(max(abs(fn-fn_rec)))

% Test exactness of SLAG transform
[f, nodes] = slag_synthesis(fn, 'P', P, 'tau', tau);
fn_rec = slag_analysis(f, 'P', P, 'tau', tau);
slag_custom_transform_error = max(max(abs(fn-fn_rec)))

nodes2 = slag_sampling(P, tau);
if (max(abs(nodes-nodes2))) ~= 0
    print('Problem with sampling scheme');
end
[f2, nodes] = slag_synthesis(fn, 'P', P, 'Nodes', nodes);
if max(max(abs(f-f2))) ~= 0
    print('Problem with transform when nodes is specified');
end

