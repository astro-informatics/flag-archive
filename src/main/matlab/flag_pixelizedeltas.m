function f = flag_pixelizedeltas( triarray, L, N, R )

sz = size(triarray);
npoints = sz(1);

[rs, thetas, phis] = flag_sampling(L, N, R);
%[rs, thetas, phis] = ndgrid(rs, thetas, phis);
%[xs, ys, zs] = ssht_s2c(thetas, phis, rs);
pixelsizes = 4*pi*( rs - [0 rs(1:N-1)] ).^2 ./ L^2 ;

f = zeros(N , L, 2*L-1);
for i= 1:npoints
    %[x, y, z] = ssht_s2c(triarray(i,2), triarray(i,3), triarray(i,1));
    %diffs = (x-xs).^2 + (y-ys).^2 + (z-zs).^2;
    %ind = find(diffs == min(min(min(diffs))));

    ind_r = find( abs(triarray(i,1)-rs) == min(abs(triarray(i,1)-rs)) );
    ind_t = find( abs(triarray(i,2)-thetas) == min(abs(triarray(i,2)-thetas)) );
    ind_p = find( abs(triarray(i,3)-phis) == min(abs(triarray(i,3)-phis)) );
    
    f(ind_r, ind_t, ind_p) = f(ind_r, ind_t, ind_p) + 1 / pixelsizes(ind_r) ;
end

end