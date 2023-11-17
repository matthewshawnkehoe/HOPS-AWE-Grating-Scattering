# test_single_eps_delta.jl
#
# DPN 6/30/21 Initial commit
# MSK 7/20/21 Added identity matrix and updated two_layer_solve and field solver
# MSK 7/26/21: Changed the size of Unm from (Nx,N+1,M+1) to (Nx,M+1,N+1)
# MSK 02/04/23: Converted from Matlab to Julia

using LinearAlgebra
using Printf
using Infiltrator

# Add the hops functions, setup functions, plot functions, sum functions, field solvers, and dno solvers
include("hops_functions.jl")
include("setup_functions.jl")
include("plot_functions.jl")
include("sum_functions.jl")
include("field_tfe_helmholtz_m_and_n.jl")
include("field_tfe_helmholtz_m_and_n_lf.jl")
include("dno_tfe_helmholtz_m_and_n.jl")
include("dno_tfe_helmholtz_m_and_n_lf.jl")
include("two_layer_solve.jl")
include("two_layer_solve_fast.jl")

# Faster Debugging
using JuliaInterpreter
using MethodAnalysis
visit(Base) do item
  isa(item, Module) && push!(JuliaInterpreter.compiled_modules, item)
  true
end

# Close all existing figures
closeall()

DoTwoLayerTest = 2
RunNumber = 1

if(RunNumber==1)
  M = 8; Nx = 16
  Eps = 1e-2; sigma = 1e-2
elseif(RunNumber==2)
  M = 6; Nx = 24
  Eps = 0.1; sigma = 0.1
elseif(RunNumber==3)
  M = 12; Nx = 32
  Eps = 0.1; sigma = 0.5
end
N = M+1; Nz = Nx

q = 1
alpha_bar = 0
d = 2*pi
c_0 = 1.0
n_u = 1.0
n_w = 1.1
# Mode = 1 (TE) or 2 (TM)
Mode = 2

A_u = 3.0
A_w = 5.0
r = 2

a = 1.0
b = 1.0
# Nz = Nx
identy = Matrix{Float64}(I, Nz+1, Nz+1)
Dz,z = cheb(Nz)

# Set up f, f_x

pp = ((2*pi/d)*[0:Nx/2-1;-Nx/2:-1]')'
xx = ((d/Nx)*(0:Nx-1)')'
f = cos.(xx)        
f_x = -sin.(xx)

err_xi_u_taylor = zeros(N+1,M+1)
err_nu_u_taylor = zeros(N+1,M+1)
err_xi_w_taylor = zeros(N+1,M+1)
err_nu_w_taylor = zeros(N+1,M+1)
err_zeta_taylor = zeros(N+1,M+1)
err_psi_taylor = zeros(N+1,M+1)
err_G_taylor = zeros(N+1,M+1)
err_J_taylor = zeros(N+1,M+1)
err_U_taylor = zeros(N+1,M+1)
err_W_taylor = zeros(N+1,M+1)
err_ubar_taylor = zeros(N+1,M+1)
err_wbar_taylor = zeros(N+1,M+1)
err_u_taylor = zeros(N+1,M+1)
err_w_taylor = zeros(N+1,M+1)

err_xi_u_pade = zeros(N+1,M+1)
err_nu_u_pade = zeros(N+1,M+1)
err_xi_w_pade = zeros(N+1,M+1)
err_nu_w_pade = zeros(N+1,M+1)
err_zeta_pade = zeros(N+1,M+1)
err_psi_pade = zeros(N+1,M+1)
err_G_pade = zeros(N+1,M+1)
err_J_pade = zeros(N+1,M+1)
err_U_pade = zeros(N+1,M+1)
err_W_pade = zeros(N+1,M+1)
err_ubar_pade = zeros(N+1,M+1)
err_wbar_pade = zeros(N+1,M+1)
err_u_pade = zeros(N+1,M+1)
err_w_pade = zeros(N+1,M+1)

err_xi_u_pade_safe = zeros(N+1,M+1)
err_nu_u_pade_safe = zeros(N+1,M+1)
err_xi_w_pade_safe = zeros(N+1,M+1)
err_nu_w_pade_safe = zeros(N+1,M+1)
err_zeta_pade_safe = zeros(N+1,M+1)
err_psi_pade_safe = zeros(N+1,M+1)
err_G_pade_safe = zeros(N+1,M+1)
err_J_pade_safe = zeros(N+1,M+1)
err_U_pade_safe = zeros(N+1,M+1)
err_W_pade_safe = zeros(N+1,M+1)
err_ubar_pade_safe = zeros(N+1,M+1)
err_wbar_pade_safe = zeros(N+1,M+1)
err_u_pade_safe = zeros(N+1,M+1)
err_w_pade_safe = zeros(N+1,M+1)

delta = sigma/(2*q+1)
  
omega_bar = q+0.5
omega = (1+delta)*omega_bar

k_u_bar = n_u*omega_bar/c_0
gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2)
_,_,alpha_bar_p,gamma_u_bar_p,eep,eem = setup_2d(Nx,d,alpha_bar,gamma_u_bar)
gamma_u_bar_r = gamma_u_bar_p[r+1]
  
k_w_bar = n_w*omega_bar/c_0
gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2)
_,_,alpha_bar_p,gamma_w_bar_p,eep,eem = setup_2d(Nx,d,alpha_bar,gamma_w_bar)
pp_r = pp[r+1]
alpha_bar_r = alpha_bar_p[r+1]
gamma_w_bar_r = gamma_w_bar_p[r+1]

k2 = alpha_bar_p[0+1]^2 + gamma_u_bar_p[0+1]^2
Tu = T_dno(alpha_bar,pp,gamma_u_bar,gamma_u_bar_p,k2,Nx,M)
  
xi_u_r_n_m, nu_u_r_n_m = setup_xi_u_nu_u_n_m(A_u,r,xx,pp,
            alpha_bar_p,gamma_u_bar_p,f,f_x,Nx,N,M)
u_n_m = field_tfe_helmholtz_m_and_n(xi_u_r_n_m,f,pp,gamma_u_bar_p,alpha_bar,
            gamma_u_bar,Dz,a,Nx,Nz,N,M,identy,alpha_bar_p)
G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m,f,pp,Dz,a,Nx,Nz,N,M)

xi_w_r_n_m, nu_w_r_n_m = setup_xi_w_nu_w_n_m(A_w,r,xx,pp,alpha_bar_p,
            gamma_w_bar_p,f,f_x,Nx,N,M)
w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_w_r_n_m,f,pp,gamma_w_bar_p,alpha_bar,
            gamma_w_bar,Dz,b,Nx,Nz,N,M,identy,alpha_bar_p)
J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m,f,pp,Dz,b,Nx,Nz,N,M)
  
if(Mode==1)
  tau2 = 1.0 # TE
else
  tau2 = (n_u/n_w)^2 # TM
end

zeta_r_n_m = xi_u_r_n_m - xi_w_r_n_m
psi_r_n_m = -nu_u_r_n_m - tau2*nu_w_r_n_m

@printf("TwoLayerSolve time: ")
if(DoTwoLayerTest==1)
  U_n_m,W_n_m,ubar_n_m,wbar_n_m = @time two_layer_solve(tau2,zeta_r_n_m,psi_r_n_m,
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,
    gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)
elseif(DoTwoLayerTest==2)
  U_n_m,W_n_m,ubar_n_m,wbar_n_m = @time two_layer_solve_fast(tau2,zeta_r_n_m,psi_r_n_m,
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,gamma_u_bar,
    gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)
else
  U_n_m = ones(Nx,M+1,N+1)
  W_n_m = U_n_m
  ubar_n_m = U_n_m
  wbar_n_m = U_n_m
end

alpha = (1+delta)*alpha_bar
gamma_u = (1+delta)*gamma_u_bar
gamma_w = (1+delta)*gamma_w_bar
alpha_r = alpha_bar_r + delta*alpha_bar
k_u = (1+delta)*k_u_bar
k_w = (1+delta)*k_w_bar
gamma_u_r = sqrt(Complex(k_u^2 - alpha_r^2))
gamma_w_r = sqrt(Complex(k_w^2 - alpha_r^2))
      
# Exact
xi_u_r = @. A_u*exp(1im*pp_r*xx)*exp(1im*gamma_u_r*Eps*f)
nu_u_r = @. (-1im*gamma_u_r + 1im*pp_r*Eps*f_x)*xi_u_r
xi_w_r = @. A_w*exp(1im*pp_r*xx)*exp(-1im*gamma_w_r*Eps*f)
nu_w_r = @. (-1im*gamma_w_r - 1im*pp_r*Eps*f_x)*xi_w_r
zeta_r = xi_u_r - xi_w_r
psi_r = -nu_u_r - tau2*nu_w_r
ubar = @. A_u*exp(1im*pp_r*xx)*exp(1im*gamma_u_r*a)
wbar = @. A_w*exp(1im*pp_r*xx)*exp(-1im*gamma_w_r*(-b))

# Exact u in transformed coordinates
ll = (0:Nz)
z_min = 0; z_max = a
tilde_z = @. cos(pi*ll/Nz)
z_prime = @. ((z_max-z_min)/2.0)*(tilde_z - 1.0) + z_max
u = zeros(Complex{Float64},Nx,Nz+1)
for j=1:Nx
  for ell=0:Nz
    z_new = @. (a-Eps*f[j])*z_prime[ell+1]/a + Eps*f[j]
    u[j,ell+1] = @. A_u*exp(1im*pp_r*xx[j])*exp(1im*gamma_u_r*z_new)
  end
end

# Exact w in transformed coordinates
ll = (0:Nz)
z_min = -b; z_max = 0
tilde_z = @. cos(pi*ll/Nz)
z_prime = @. ((z_max-z_min)/2.0)*(tilde_z - 1.0) + z_max
w = zeros(Complex{Float64},Nx,Nz+1)
for j=1:Nx
  for ell=0:Nz
    z_new = @. (b+Eps*f[j])*z_prime[ell+1]/b + Eps*f[j]
    w[j,ell+1] = @. A_w*exp(1im*pp_r*xx[j])*exp(-1im*gamma_w_r*z_new)
  end
end

@printf("Summation time: ")
@time for n=0:N
  for m=0:M
    # Taylor
    xi_u_approx = fcn_sum(1,xi_u_r_n_m,Eps,delta,Nx,n,m)
    nu_u_approx = fcn_sum(1,nu_u_r_n_m,Eps,delta,Nx,n,m)
    xi_w_approx = fcn_sum(1,xi_w_r_n_m,Eps,delta,Nx,n,m)
    nu_w_approx = fcn_sum(1,nu_w_r_n_m,Eps,delta,Nx,n,m)
    zeta_approx = fcn_sum(1,zeta_r_n_m,Eps,delta,Nx,n,m)
    psi_approx = fcn_sum(1,psi_r_n_m,Eps,delta,Nx,n,m)
    G_approx = fcn_sum(1,G_n_m,Eps,delta,Nx,n,m)
    J_approx = fcn_sum(1,J_n_m,Eps,delta,Nx,n,m)
    U_approx = fcn_sum(1,U_n_m,Eps,delta,Nx,n,m)
    W_approx = fcn_sum(1,W_n_m,Eps,delta,Nx,n,m)
    ubar_approx = fcn_sum(1,ubar_n_m,Eps,delta,Nx,n,m)
    wbar_approx = fcn_sum(1,wbar_n_m,Eps,delta,Nx,n,m)
    u_approx = vol_fcn_sum(1,u_n_m,Eps,delta,Nx,Nz,n,m)
    w_approx = vol_fcn_sum(1,w_n_m,Eps,delta,Nx,Nz,n,m)
    
    err_xi_u_taylor[n+1,m+1] = norm(xi_u_r-xi_u_approx,Inf)
    err_nu_u_taylor[n+1,m+1] = norm(nu_u_r-nu_u_approx,Inf)
    err_xi_w_taylor[n+1,m+1] = norm(xi_w_r-xi_w_approx,Inf)
    err_nu_w_taylor[n+1,m+1] = norm(nu_w_r-nu_w_approx,Inf)
    err_zeta_taylor[n+1,m+1] = norm(zeta_r-zeta_approx,Inf)
    err_psi_taylor[n+1,m+1] = norm(psi_r-psi_approx,Inf)
    err_G_taylor[n+1,m+1] = norm(nu_u_r-G_approx,Inf)
    err_J_taylor[n+1,m+1] = norm(nu_w_r-J_approx,Inf)
    err_U_taylor[n+1,m+1] = norm(xi_u_r-U_approx,Inf)
    err_W_taylor[n+1,m+1] = norm(xi_w_r-W_approx,Inf)
    err_ubar_taylor[n+1,m+1] = norm(ubar-ubar_approx,Inf)
    err_wbar_taylor[n+1,m+1] = norm(wbar-wbar_approx,Inf)
    err_u_taylor[n+1,m+1] = maximum(abs.(u-u_approx))
    err_w_taylor[n+1,m+1] = maximum(abs.(u-u_approx))
  
    # Pade  
    xi_u_approx = fcn_sum(2,xi_u_r_n_m,Eps,delta,Nx,n,m)
    nu_u_approx = fcn_sum(2,nu_u_r_n_m,Eps,delta,Nx,n,m)
    xi_w_approx = fcn_sum(2,xi_w_r_n_m,Eps,delta,Nx,n,m)
    nu_w_approx = fcn_sum(2,nu_w_r_n_m,Eps,delta,Nx,n,m)
    zeta_approx = fcn_sum(2,zeta_r_n_m,Eps,delta,Nx,n,m)
    psi_approx = fcn_sum(2,psi_r_n_m,Eps,delta,Nx,n,m)
    G_approx = fcn_sum(2,G_n_m,Eps,delta,Nx,n,m)
    J_approx = fcn_sum(2,J_n_m,Eps,delta,Nx,n,m)
    U_approx = fcn_sum(2,U_n_m,Eps,delta,Nx,n,m)
    W_approx = fcn_sum(2,W_n_m,Eps,delta,Nx,n,m)
    ubar_approx = fcn_sum(2,ubar_n_m,Eps,delta,Nx,n,m)
    wbar_approx = fcn_sum(2,wbar_n_m,Eps,delta,Nx,n,m)
    u_approx = vol_fcn_sum(2,u_n_m,Eps,delta,Nx,Nz,n,m)
    w_approx = vol_fcn_sum(2,w_n_m,Eps,delta,Nx,Nz,n,m)
    
    err_xi_u_pade[n+1,m+1] = norm(xi_u_r-xi_u_approx,Inf)
    err_nu_u_pade[n+1,m+1] = norm(nu_u_r-nu_u_approx,Inf)
    err_xi_w_pade[n+1,m+1] = norm(xi_w_r-xi_w_approx,Inf)
    err_nu_w_pade[n+1,m+1] = norm(nu_w_r-nu_w_approx,Inf)
    err_zeta_pade[n+1,m+1] = norm(zeta_r-zeta_approx,Inf)
    err_psi_pade[n+1,m+1] = norm(psi_r-psi_approx,Inf)
    err_G_pade[n+1,m+1] = norm(nu_u_r-G_approx,Inf)
    err_J_pade[n+1,m+1] = norm(nu_w_r-J_approx,Inf)
    err_U_pade[n+1,m+1] = norm(xi_u_r-U_approx,Inf)
    err_W_pade[n+1,m+1] = norm(xi_w_r-W_approx,Inf)
    err_ubar_pade[n+1,m+1] = norm(ubar-ubar_approx,Inf)
    err_wbar_pade[n+1,m+1] = norm(wbar-wbar_approx,Inf)
    err_u_pade[n+1,m+1] = maximum(abs.(u-u_approx))
    err_w_pade[n+1,m+1] = maximum(abs.(u-u_approx))
    
    # Pade (safe)
    xi_u_approx = fcn_sum(3,xi_u_r_n_m,Eps,delta,Nx,n,m)
    nu_u_approx = fcn_sum(3,nu_u_r_n_m,Eps,delta,Nx,n,m)
    xi_w_approx = fcn_sum(3,xi_w_r_n_m,Eps,delta,Nx,n,m)
    nu_w_approx = fcn_sum(3,nu_w_r_n_m,Eps,delta,Nx,n,m)
    zeta_approx = fcn_sum(3,zeta_r_n_m,Eps,delta,Nx,n,m)
    psi_approx = fcn_sum(3,psi_r_n_m,Eps,delta,Nx,n,m)
    G_approx = fcn_sum(3,G_n_m,Eps,delta,Nx,n,m)
    J_approx = fcn_sum(3,J_n_m,Eps,delta,Nx,n,m)
    U_approx = fcn_sum(3,U_n_m,Eps,delta,Nx,n,m)
    W_approx = fcn_sum(3,W_n_m,Eps,delta,Nx,n,m)
    ubar_approx = fcn_sum(3,ubar_n_m,Eps,delta,Nx,n,m)
    wbar_approx = fcn_sum(3,wbar_n_m,Eps,delta,Nx,n,m)
    u_approx = vol_fcn_sum(3,u_n_m,Eps,delta,Nx,Nz,n,m)
    w_approx = vol_fcn_sum(3,w_n_m,Eps,delta,Nx,Nz,n,m)
    
    err_xi_u_pade_safe[n+1,m+1] = norm(xi_u_r-xi_u_approx,Inf)
    err_nu_u_pade_safe[n+1,m+1] = norm(nu_u_r-nu_u_approx,Inf)
    err_xi_w_pade_safe[n+1,m+1] = norm(xi_w_r-xi_w_approx,Inf)
    err_nu_w_pade_safe[n+1,m+1] = norm(nu_w_r-nu_w_approx,Inf)
    err_zeta_pade_safe[n+1,m+1] = norm(zeta_r-zeta_approx,Inf)
    err_psi_pade_safe[n+1,m+1] = norm(psi_r-psi_approx,Inf)
    err_G_pade_safe[n+1,m+1] = norm(nu_u_r-G_approx,Inf)
    err_J_pade_safe[n+1,m+1] = norm(nu_w_r-J_approx,Inf)
    err_U_pade_safe[n+1,m+1] = norm(xi_u_r-U_approx,Inf)
    err_W_pade_safe[n+1,m+1] = norm(xi_w_r-W_approx,Inf)
    err_ubar_pade_safe[n+1,m+1] = norm(ubar-ubar_approx,Inf)
    err_wbar_pade_safe[n+1,m+1] = norm(wbar-wbar_approx,Inf)
    err_u_pade_safe[n+1,m+1] = maximum(abs.(u-u_approx))
    err_w_pade_safe[n+1,m+1] = maximum(abs.(u-u_approx))
  end
end

# Make error plots
t1,t2,t3,t4,t5,t6 = plot_errors("Taylor",N,M,err_xi_u_taylor,err_nu_u_taylor,err_zeta_taylor,
        err_xi_w_taylor,err_nu_w_taylor,err_psi_taylor)
p1,p2,p3,p4,p5,p6 = plot_errors("Pade",N,M,err_xi_u_pade,err_nu_u_pade,err_zeta_pade,
        err_xi_w_pade,err_nu_w_pade,err_psi_pade)
ps1,ps2,ps3,ps4,ps5,ps6 = plot_errors("Pade Safe",N,M,err_xi_u_pade_safe,err_nu_u_pade_safe,err_zeta_pade_safe,
        err_xi_w_pade_safe,err_nu_w_pade_safe,err_psi_pade_safe)
t7,t8,t9,t10,t11,t12 = plot_errors("Taylor",N,M,err_G_taylor,err_U_taylor,err_ubar_taylor,
        err_J_taylor,err_W_taylor,err_wbar_taylor)
p7,p8,p9,p10,p11,p12= plot_errors("Pade",N,M,err_G_pade,err_U_pade,err_ubar_pade,
        err_J_pade,err_W_pade,err_wbar_pade)
ps7,ps8,ps9,ps10,ps11,ps12 = plot_errors("Pade Safe",N,M,err_G_pade_safe,err_U_pade_safe,err_ubar_pade_safe,
        err_J_pade_safe,err_W_pade_safe,err_wbar_pade_safe)
f1,f2,f3,f4,f5,f6 = plot_errors("Fields",N,M,err_u_taylor,err_u_pade,err_u_pade_safe,
        err_w_taylor,err_w_pade,err_w_pade_safe)

# Display error plots on external windows      
fig1 = plot(t1,t2,t3,t4,t5,t6,layout=(2,3),margin=2mm,dpi=200)
fig2 = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),margin=2mm,dpi=200)
fig3 = plot(ps1,ps2,ps3,ps4,ps5,ps6,layout=(2,3),margin=2mm,dpi=200)
fig4 = plot(t7,t8,t9,t10,t11,t12,layout=(2,3),margin=2mm,dpi=200)
fig5 = plot(p7,p8,p9,p10,p11,p12,layout=(2,3),margin=2mm,dpi=200)
fig6 = plot(ps7,ps8,ps9,ps10,ps11,ps12,layout=(2,3),margin=2mm,dpi=200)
fig7 = plot(f1,f2,f3,f4,f5,f6,layout=(2,3),margin=2mm,dpi=200)
display(fig1)
display(fig2)
display(fig3)
display(fig4)
display(fig5)
display(fig6)
display(fig7)