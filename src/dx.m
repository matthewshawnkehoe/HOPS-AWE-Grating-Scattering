function [u_x] = dx(u,p)

u_x = ifft((1i*p).*fft(u));

return;