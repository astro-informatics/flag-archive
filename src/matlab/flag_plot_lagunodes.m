function flag_plot_lagunodes( P, R )

figure('Position',[200 200 800 800])
hold on
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
nodes = R;
for i=2:P
    nodes = slag_sampling(i, R);
    y = i * nodes ./ nodes ;
    plot( nodes, y, '--', 'LineWidth',5, 'Color', 'black' )
    plot( nodes, y, 'o','MarkerSize',15, 'LineWidth',5, 'Color', 'black' )
end
axis([0 R 1 P+1])
set(gca,'FontSize',42);
set(gca,'LineWidth',5);

end