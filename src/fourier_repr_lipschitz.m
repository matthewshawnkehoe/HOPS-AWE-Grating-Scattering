function [f,f_x] = fourier_repr_lipschitz(P,x)

f = 0;
for k = 1:ceil(P/2)
    f = f + 8/(pi^2*(2*k-1)^2) * cos((2*k-1)*x);
end

f_x = 0;
for k = 1:ceil(P/2)
    f_x = f_x + (-8*(2*k-1))/(pi^2*(2*k-1)^2) * sin((2*k-1)*x);
end

figure(3);
hold on;
subplot(1,2,1)
plot(x,f);
title('Lipschitz profile for $f$','interpreter','latex','FontSize',16);
subplot(1,2,2)
plot(x,f_x);
title('Lipschitz profile for $f_x$','interpreter','latex','FontSize',16);

return