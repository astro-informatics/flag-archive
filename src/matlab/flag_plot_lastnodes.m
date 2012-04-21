% flag_plot_lastnodes
% Show a good approximation for the last node of laguerre scheme

nmax = 400;
alpha = 2;
nrange = 20:nmax;
lastnodes = zeros(1, max(size(nrange)));
for N = nrange;
   [nodes, weights] = slag_gausslaguerre_quadrature(N+1, alpha);
   lastnodes(N-nrange(1)+1) = nodes(N+1);
end


plot(nrange, lastnodes, '--o', 'linewidth', 2)

hold on
law = 3.95*nrange + 10;
plot(nrange, law, 'r', 'linewidth', 2)