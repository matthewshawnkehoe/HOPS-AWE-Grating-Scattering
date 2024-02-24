% Constants
a = 1.0;
alpha = 0;
gamma = 1.5;
Nx = 32;
Nz = 32;
N = 20;
M = 20;

% N-Dimensional Test Matrices
p = rand(Nx, 1);
real_part = rand(Nx, Nz + 1, N + 1, M + 1);
imaginary_part = rand(Nx, Nz + 1, N + 1, M + 1);
unm = complex(real_part, imaginary_part);
A1_xx = rand(Nx, Nz + 1);
A1_xz = rand(Nx, Nz + 1);
A1_zx = rand(Nx, Nz + 1);
A2_xx = rand(Nx, Nz + 1);
A2_xz = rand(Nx, Nz + 1);
A2_zx = rand(Nx, Nz + 1);
A2_zz = rand(Nx, Nz + 1);
B1_x = rand(Nx, Nz + 1);
B2_x = rand(Nx, Nz + 1);
B2_z = rand(Nx, Nz + 1);
Dz = rand(Nx + 1, Nz + 1);
S1 = rand(Nx, Nz + 1);
S2 = rand(Nx, Nz + 1);

% Create grids for n and m
[nGrid, mGrid] = ndgrid(0:N, 0:M);

% Compute u_x_all and u_z_all for all combinations of n and m
[u_x_all, u_z_all] = deal(zeros(Nx, Nz + 1, N + 1, M + 1));
[u_x_all(:, :, nGrid + 1, mGrid + 1), u_z_all(:, :, nGrid + 1, mGrid + 1)] = deal(dx(unm(:, :, mGrid + 1, nGrid + 1), p), dz(unm(:, :, mGrid + 1, nGrid + 1), Dz, a));

% Compute Fnm using vectorized operations
Fnm = zeros(Nx, Nz + 1);
if any(nGrid(:) >= 1)
    temp = A1_xx .* u_x_all(:, :, :, nGrid + 1) ...
         - dx(A1_zx .* u_x_all(:, :, :, nGrid + 1), p) ...
         - dz(A1_xx .* u_z_all(:, :, :, nGrid + 1), Dz, a) ...
         - B1_x .* u_x_all(:, :, :, nGrid + 1) ...
         - dx(A1_xz .* u_z_all(:, :, :, nGrid + 1), p) ...
         - 2 * 1i * alpha * S1 .* u_x_all(:, :, :, nGrid + 1) ...
         - gamma^2 * S1 .* unm(:, :, :, nGrid + 1);
    Fnm = Fnm - sum(temp, 4);
end

if any(mGrid(:) >= 1)
    temp = 2 * 1i * alpha * u_x_all(:, :, :, mGrid + 1) ...
         - 2 * gamma^2 * unm(:, :, :, mGrid + 1);
    Fnm = Fnm - sum(temp, 4);
end

if any(nGrid(:) >= 1 & mGrid(:) >= 1)
    temp = 2 * 1i * alpha * S2 .* u_x_all(:, :, :, mGrid + 1) ...
         - 2 * gamma^2 * S2 .* unm(:, :, :, mGrid + 1);
    Fnm = Fnm - sum(temp, 4);
end

if any(nGrid(:) >= 2)
    temp = A2_xx .* u_x_all(:, :, :, nGrid + 1) ...
         - dx(A2_zx .* u_x_all(:, :, :, nGrid + 1), p) ...
         - dz(A2_xx .* u_z_all(:, :, :, nGrid + 1), Dz, a) ...
         - B2_x .* u_x_all(:, :, :, nGrid + 1) ...
         - dx(A2_xz .* u_z_all(:, :, :, nGrid + 1), p) ...
         - gamma^2 * S2 .* unm(:, :, :, nGrid - 1 + 1);
    Fnm = Fnm - sum(temp, 4);
end

if any(mGrid(:) >= 2)
    temp = gamma^2 * unm(:, :, :, mGrid - 1 + 1);
    Fnm = Fnm - sum(temp, 4);
end

if any(nGrid(:) >= 1 & mGrid(:) >= 2)
    temp = 2 * 1i * alpha * S2 .* u_x_all(:, :, :, mGrid - 1 + 1) ...
         - 2 * gamma^2 * S2 .* unm(:, :, :, mGrid - 1 + 1);
    Fnm = Fnm - sum(temp, 4);
end

% Functions
function u_x = dx(u, p)
u_x = ifft((1i * p) .* fft(u), [], 1);
end

function u_z = dz(u, Dz, b)
    % Compute u_z using matrix multiplication and appropriate reshaping
    u_z = reshape(((2.0 / b) * reshape(Dz, [], size(Dz, 2)) * reshape(u, size(u, 1), [])), size(u));
end
