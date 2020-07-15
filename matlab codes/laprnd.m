function y  = laprnd(mu, sigma, m, n)
b = sigma/sqrt(2);
u = rand(m, n)-0.5;
y = mu - b*sign(u).* log(1- 2*abs(u));
end 