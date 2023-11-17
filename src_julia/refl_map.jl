# refl_map.m

using Printf
using Profile
using Base.Threads
using BenchmarkTools
#using PlotlyJS

# Add the hops functions, setup functions, plot functions, field solvers, dno solvers, two layer solve, and enery defect
include("hops_functions.jl")
include("setup_functions.jl")
include("plot_functions.jl")
include("field_tfe_helmholtz_m_and_n.jl")
include("field_tfe_helmholtz_m_and_n_lf.jl")
include("dno_tfe_helmholtz_m_and_n.jl")
include("dno_tfe_helmholtz_m_and_n_lf.jl")
include("two_layer_solve_fast.jl")
include("energy_defect.jl")

# Faster Debugging
using JuliaInterpreter
using MethodAnalysis
visit(Base) do item
  isa(item, Module) && push!(JuliaInterpreter.compiled_modules, item)
  true
end

# Close all existing figures
closeall()

PlotLambda = 1
PlotRelative = 1

RunNumber = 100

if(RunNumber==1)
  M = 4; Nx = 16
  Eps_Max = 1e-2; sigma = 1e-2
elseif(RunNumber==2)
  M = 6; Nx = 24
  Eps_Max = 0.1; sigma = 0.1
elseif(RunNumber==3)
  M = 8; Nx = 32
  Eps_Max = 0.1; sigma = 0.5
elseif(RunNumber==100)
  # HOPS/AWE paper
  M = 15; Nx = 32
  Eps_Max = 0.2; sigma = 0.99
end
N = M; Nz = 32

alpha_bar = 0
d = 2*pi
c_0 = 1
n_u = 1.0
n_w = 0.05 + 2.275*1im   #2.3782, Carbon
# n_w = 1.1
# Mode = 1 (TE) or 2 (TM)
Mode = 2
Taylor = false

N_delta = 100
N_Eps = 100
# qq = (1:2)'
# qq = (1)
qq = (1:6)'

a = 1.0; b = 1.0
# Nz = Nx
identy = Matrix{Float64}(I, Nz+1, Nz+1)
Dz,z = cheb(Nz)

Eps = LinRange(0,Eps_Max,N_Eps)

# Set up f, f_x

pp = ((2*pi/d)*[0:Nx/2-1;-Nx/2:-1]')'
xx = ((d/Nx)*(0:Nx-1)')'
#f =  (1/4)*sin.(4*pi*xx/d)
#f_x = (pi/d)*cos.(4*pi*xx/d)
#f = (1/2)*cos.(2*pi*xx/d)
#f_x = -(2*pi/d)*sin.(2*pi*xx/d)
f = cos.(4*xx)
f_x = -4*sin.(4*xx)
# f = cos.(xx)
# f_x = -sin.(xx)

# Sawtooth/Lipschitz
# P = 40
# f,f_x = fourier_repr_lipschitz(P,xx)

# Rough Profile
# P = 40
# f,f_x = fourier_repr_rough(P,xx)

R_taylor = zeros(N_Eps,N_delta)
R_pade = zeros(N_Eps,N_delta)
R_pade_safe = zeros(N_Eps,N_delta)

# Create plot figures
if(PlotLambda==0)
  c_ed = contourf(title=L"D",xaxis=L"\omega",yaxis=L"\varepsilon",color=:hot)
  c_rm = contourf(title=L"R",xaxis=L"\omega",yaxis=L"\varepsilon",color=:hot)
else
  c_ed = contourf(title=L"D",xaxis=L"\lambda",yaxis=L"\varepsilon",color=:hot)
  c_rm = contourf(title=L"R",xaxis=L"\lambda",yaxis=L"\varepsilon",color=:hot)
end

# Loop over q
# TODO - Enable multi-threading through adding @threads below?
for s in eachindex(qq)
  q = qq[s]
  if(N_delta==1)
    delta = [0]
  else
    delta = LinRange(-sigma/(2*q+1),sigma/(2*q+1),N_delta)
  end
  
  omega_bar = c_0*(2*pi/d)*(q + 0.5)
  omega = @. (1+delta)*omega_bar
  lambda = @. (2*pi*c_0/omega)

  k_u_bar = n_u*omega_bar/c_0
  gamma_u_bar = sqrt(k_u_bar^2 - alpha_bar^2)
  _,_,alpha_bar_p,gamma_u_bar_p,eep,eem = setup_2d(Nx,d,alpha_bar,gamma_u_bar)

  k_w_bar = n_w*omega_bar/c_0
  gamma_w_bar = sqrt(k_w_bar^2 - alpha_bar^2)
  _,_,alpha_bar_p,gamma_w_bar_p,eep,eem = setup_2d(Nx,d,alpha_bar,gamma_w_bar)

  zeta_n_m,psi_n_m = setup_zeta_psi_n_m(xx,pp,alpha_bar,gamma_u_bar,f,f_x,Nx,N,M)
  
  if(Mode==1)
    tau2 = 1.0 # TE
  else
    tau2 = (n_u/n_w)^2 # TM
  end

  @printf("TwoLayerSolve time: ")
  U_n_m,W_n_m,ubar_n_m,wbar_n_m = @time two_layer_solve_fast(tau2,zeta_n_m,psi_n_m,
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,
    gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)

  # Correct order of dimensions
  ubar_n_m_upd = permutedims(ubar_n_m,[3 2 1])
  wbar_n_m_upd = permutedims(wbar_n_m,[3 2 1])
  
  ee_flat,ru_flat,rl_flat = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,
      d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,
      Nx,0,0,N_Eps,N_delta,1)
  
  # Don't compute pade_sum unless Taylor is false
  if Taylor == true
    ee_taylor,ru_taylor,rl_taylor = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,
        d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,
        Nx,N,M,N_Eps,N_delta,1)
  else
    ee_pade,ru_pade,rl_pade = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,
        d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,
        Nx,N,M,N_Eps,N_delta,2)
  end
  
  # We aren't currently testing pade_safe
  # ee_pade_safe,ru_pade_safe,rl_pade_safe = energy_defect(tau2,ubar_n_m_upd,wbar_n_m_upd,
  #     d,alpha_bar,gamma_u_bar,gamma_w_bar,Eps,delta,
  #     Nx,N,M,N_Eps,N_delta,3)

  # Plot the energy defect (log10) and reflectivity map
  if Taylor == true
    plot_energy_defect_taylor(c_ed,PlotLambda,omega,lambda,Eps,ee_taylor)
    plot_refl_map_taylor(c_rm,PlotRelative,omega,lambda,Eps,ru_taylor,ru_flat)
  else
    plot_energy_defect_pade(c_ed,PlotLambda,omega,lambda,Eps,ee_pade)
    plot_refl_map_pade(c_rm,PlotRelative,omega,lambda,Eps,ru_pade,ru_flat)
  end

  # Display the contour plots
  display(c_ed)
  display(c_rm)
end
