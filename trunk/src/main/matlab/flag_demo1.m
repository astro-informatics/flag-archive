function flag_demo1

P = 10    % Number of points / band-limit for the radial sampling
npoints = 64   % Oversampling the radial basis functions
R = 1.0   % Radial size of the ball
L = 10    % Angular band-limit

flag_plot_lagunodes( P, R )

flag_plot_laguerre_kernels(P, R, npoints)

flag_plot_sampling(L, P, R)

end