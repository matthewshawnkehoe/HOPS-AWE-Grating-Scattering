using FFTW

# Add fcn_sum_fast from sum_functions.jl
include("sum_functions.jl")

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
      #println(ru[j,ell])
      rl[j,ell] = tau2*sum(PropMode_w.*C)
      #println(rl[j,ell])
    end
  end

  ee = @. 1.0 - ru - rl

  return ee,ru,rl
end