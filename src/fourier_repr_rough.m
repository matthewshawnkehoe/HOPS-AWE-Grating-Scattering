function [f,f_x] = fourier_repr_rough(P,x)

f = 0;
for k = 1:P
    f = f + 96*(2*k^2*pi^2 - 21)/(125*k^8) * cos(k*x);
end

f_x = 0;
for k = 1:P
    f_x = f_x + (-96*k*(2*k^2*pi^2 - 21))/(125*k^8) * sin(k*x);
end

figure(4);
hold on;
subplot(1,2,1)
plot(x,f);
title('Rough profile for $f$','interpreter','latex','FontSize',16);
subplot(1,2,2)
plot(x,f_x);
title('Rough profile for $f_x$','interpreter','latex','FontSize',16);

return