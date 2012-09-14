function flag_plot_laguerre_kernels(P, R, npoints)

h = R/(npoints);
nodes = (0:h:R);
sampling = slag_sampling(P, R);
            
nmax = P;%10;
colors = rand(P,3)*0.9;
% 
% figure('Position',[1 1 500 500])
% hold on
% for n = 1:nmax
%     fn = zeros(1,P);
%     fn(n) = 1.0;
%     [f, nodes] = slag_synthesis(fn, P, 'Podes', nodes);
%     p = plot(nodes, f.*nodes/sqrt(tau));
%     set(p,'Color',colors(n,:),'LineWidth',1.1);
% end
% plot(sampling, 0*sampling,'d','LineWidth',1.1,'MarkerSize',10,'MarkerEdgeColor','k',...
%                 'MarkerFaceColor',[.49 1 .63])
% axis([0 R -1 1])
% xlabel('Radius')
% title('Decaying Laguerre polynomials (weighted with exp-r)')

figure('Position',[200 200 600 1000])

subplot(3,1,1)
hold on
for n = 1:nmax
    fn = zeros(1,P);
    fn(n) = 1.0;
    [f, nodes] = slag_synthesis(fn, 'P', P, 'Nodes', nodes);
    p = plot(nodes, f);
    set(p,'Color',colors(n,:),'LineWidth',2);
end
p=plot(sampling, 0*sampling,'.');
set(p,'Color','black','LineWidth',2);
plot(sampling, 0*sampling,'o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k')
axis([0 R -0.3 2.6])
xlabel('x1','FontSize',10)
ylabel('y1','FontSize',10)
set(gca,'FontSize',10);
set(gca,'LineWidth',2);

subplot(3,1,2)
hold on
for n = 1:nmax
    fn = zeros(1,P);
    fn(n) = 1.0;
    [f, nodes] = slag_synthesis(fn, 'P', P, 'Nodes', nodes);
    p = plot(nodes, f);
    set(p,'Color',colors(n,:),'LineWidth',2);
end
p=plot(sampling, 0*sampling,'.');
set(p,'Color','black','LineWidth',2);
plot(sampling, 0*sampling,'o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k')
axis([0 R -0.2 0.2])
xlabel('x2','FontSize',10)
ylabel('y2','FontSize',10)
set(gca,'FontSize',10);
set(gca,'LineWidth',2);

subplot(3,1,3)
hold on
for n = 1:nmax
    fn = zeros(1,P);
    fn(n) = 1.0;
    [f, nodes] = slag_synthesis(fn, 'P', P, 'Nodes', nodes);
    p = plot(nodes, f.*nodes);
    set(p,'Color',colors(n,:),'LineWidth',2);
end
p=plot(sampling, 0*sampling,'.');
set(p,'Color','black','LineWidth',2);
plot(sampling, 0*sampling,'o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k')
axis([0 R -0.029 0.029])
xlabel('x3','FontSize',10)
ylabel('y3','FontSize',10)
set(gca,'FontSize',10);
set(gca,'LineWidth',2);

end