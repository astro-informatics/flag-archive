function flag_bessel_demo2(P, ell)

coefs = zeros(P,P);
for p = 0:P-1
    cjp = 2 * factorial(ell+2) * (p+1)*(p+2);
    coefs(p+1,1) = cjp;
    for j = 1:p
        cjp = (p-j+1)/(j*(j+2)) * (j+ell+2) * 2 * cjp ;% * (-1) ; 
        coefs(p+1,j+1) = cjp;
    end
end


log(coefs)
plot(log(coefs))

end