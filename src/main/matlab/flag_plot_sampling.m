function flag_plot_sampling(L, P, R)

[rs, thetas, phis] = flag_sampling(L, P, R);

x = [];
y = [];
z = [];
s = [];
for i = 1:P
    for j = 1:L
        for k = 1:2*L-1
            x = [x rs(i)*cos(phis(k))*sin(thetas(j))];
            y = [y rs(i)*sin(phis(k))*sin(thetas(j))];
            z = [z rs(i)*cos(thetas(j))];
            s = [s i];
        end
    end   
end
c = numel(x):-1:1;

figure('Position',[1 1 600 600])
h = scatter3(x,y,z,0.5*s,c,'filled');
set(gca, 'visible', 'off'); 
view(60,30);
zoom(1.5);
end