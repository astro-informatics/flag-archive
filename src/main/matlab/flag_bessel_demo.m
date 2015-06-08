function flag_bessel_demo(N, Nk)

R = 200;
ell = 0; % Careful! only true for ell=0!
kmin=0.05;
kmax=0.5;
kh=(kmax-kmin)/Nk;
kvalues = kmin:kh:kmax;

tau = flag_get_tau(N, R);
nodes = slag_sampling(N, R);
%nodes = nodes/tau
%tau=1;

f = nodes.^(-2);
fn = slag_analysis(f);
frec = slag_synthesis(fn);

flk_true = sqrt(pi/2) ./ kvalues;


%syms r;
Rmax = 5000;

flk_num = zeros(size(flk_true));
flk = zeros(size(flk_true));
for k = 1:numel(kvalues)
    ktilde = tau * kvalues(k);
    z = -4*ktilde^2;
    jlpk = zeros(1,N);
    mujlk = zeros(1,N+2);
    mujlk_bis = zeros(1,N+2);
    for j = 0:numel(fn)+1
       hypervalue = hypergeom([(j+ell+1)/2.0, (j+ell)/2.0+1], ell+1.5, -4*ktilde^2) / gamma(ell+1.5);
       %[j, (j+ell+1)/2.0, (j+ell)/2.0+1, ell+1.5, -4*ktilde^2, hypervalue]
       %fbis = exp(-r/(2*tau)) * r^j * sin(kvalues(k)*r) / (kvalues(k)*r);
       %mujlk_bis(1,j+1) =  vpa(int(fbis, r, 0, Rmax)) * tau^(-j+2.0);  
       mujlk(1,j+1) = sqrt(pi) * ktilde^ell * tau^3 * hypervalue;
    end
    %kvalues(k)
    %tau
    %mujlk
    for p = 0:numel(fn)-1
        cjp = 2 * factorial(ell+2) * (p+1)*(p+2);
        jlpk(1,p+1) = ((p+1)*(p+2))^(-0.5) * cjp * mujlk(3);
        for j = 1:p
            cjp = -(p-j+1)/(j*(j+2)) * (j+ell+2) * 2 * cjp; 
            jlpk(1,p+1) = jlpk(1,p+1) + ((p+1)*(p+2))^(-0.5) * cjp * mujlk(j+3);
        end
    end
    jlpk_bis = zeros(1,N);
    for p = 0:numel(fn)-1
        kp = @(r) exp(-r./(2.*tau)) .* r.^2 .* mfun('L', p, 2, r./tau) .* sin(kvalues(k).*r) ./ (kvalues(k).*r);
        jlpk_bis(1,p+1) = ((p+1)*(p+2))^(-0.5) .* quad(kp,0,Rmax);
    end
    jlpk
    jlpk_bis
    flk(k) = sqrt(2/pi) * sum( fn .* jlpk );
    %fbis = @(r) sin(kvalues(k).*r) ./ (kvalues(k).*r);
    %flk_num(k) = sqrt(2/pi) * quad(fbis, 0, Rmax);
end

flk_true
flk
%flk_num
flk_true./flk

end