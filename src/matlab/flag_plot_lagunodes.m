function flag_plot_lagunodes( P, R )

figure('Position',[200 200 900 900])
hold on
title('Rescaled Laguerre sampling','FontSize',42)
xlabel('Radius','FontSize',42)
ylabel('Number of samples','FontSize',42)
nodes = R;
for i=2:P
    nodes = slag_sampling(i, R);
    y = i * nodes ./ nodes ;
    plot( nodes, y, '--', 'LineWidth',5, 'Color', 'black' )
    plot( nodes, y, 'o','MarkerSize',20, 'LineWidth',5, 'Color', 'black' )
end
axis([0 R 1 P+1])
set(gca,'FontSize',42);
set(gca,'LineWidth',5);

end