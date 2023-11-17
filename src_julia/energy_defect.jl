using FFTW

# Add fcn_sum_fast from sum_functions.jl
include("sum_functions.jl")

"""
Compute the energy defect (D) and inspect the conservation of energy.
### Input
- `tau2` --  a numerical constant representing TE or TM mode
- `ubar_n_m` -- a tensor representing approximate solution in upper layer
- `wbar_n_m` --  tensor representing approximate solution in lower layer
- `d` -- the periodicity of the grating
- `alpha_bar` -- a numerical constant
- `gamma_u_bar` -- a numerical constant in the upper field
- `gamma_w_bar` -- a numerical constant in the lower field
- `Eps` -- the physical error in the surface deformation
- `delta` -- the numerical error in the discretization of the frequency perbutation
- `N` -- the maximum number of Taylor orders for the interfacial perturbation
- `M` -- the maximum number of Taylor orders for the frequency perturbation
- `N_Eps` -- the maximum number of epsilon orders for the interfacial perturbation
- `N_Delta` -- the maximum number of delta orders for the frequency perturbation
- `SumType` -- boolean to control Taylor or Pade summation 
### Output
- `ee` -- an estimate of the total energy defect (which should approach 1)
- `ru` -- proportion of scattered energy in the upper field
- `rl` -- proportion of scattered energy in the lower field
"""
function energy_defect(tau2,ubar_n_m,wbar_n_m, d,alpha_bar,gamma_u_bar,
    gamma_w_bar,Eps,delta,Nx,N,M,N_Eps,N_delta,SumType)

  ru = zeros(N_Eps,N_delta)
  rl = zeros(N_Eps,N_delta)

  k_u_bar = sqrt(alpha_bar^2 + gamma_u_bar^2)
  k_w_bar = sqrt(alpha_bar^2 + gamma_w_bar^2)

  PropMode_u = zeros(Nx,1)
  PropMode_w = zeros(Nx,1)

  for ell=1:N_delta
    for j=1:N_Eps
      alpha = (1+delta[ell])*alpha_bar
      gamma_u = (1+delta[ell])*gamma_u_bar
      gamma_w = (1+delta[ell])*gamma_w_bar
      k_u = (1+delta[ell])*k_u_bar
      k_w = (1+delta[ell])*k_w_bar
    
      xx,pp,alpha_p,gamma_u_p,eep,eem = setup_2d(Nx,d,alpha,gamma_u)
      xx,pp,alpha_p,gamma_w_p,eep,eem = setup_2d(Nx,d,alpha,gamma_w)
    
      gamma_i = gamma_u_p[0+1]
      
      ubar = fcn_sum_fast(SumType,ubar_n_m,Eps[j],delta[ell],Nx,N,M)
      wbar = fcn_sum_fast(SumType,wbar_n_m,Eps[j],delta[ell],Nx,N,M)
  
      ubarhat = fft(ubar,(1,))/Nx
      wbarhat = fft(wbar,(1,))/Nx
      
      alpha_p_2 = alpha_p.^2
      k_u_2 = k_u^2
      k_w_2 = k_w^2
      
      B = @. gamma_u_p/gamma_i*abs(ubarhat)^2
      C = @. gamma_w_p/gamma_i*abs(wbarhat)^2
      
      @. PropMode_u[alpha_p_2 < k_u_2] = 1
      @. PropMode_w[alpha_p_2 < real(k_w_2)] = 1

      ru[j,ell] = sum(PropMode_u.*B)
      rl[j,ell] = tau2*sum(PropMode_w.*C)
    end
  end

  ee = @. 1.0 - ru - rl

  return ee,ru,rl
end