# Add the hops functions, setup functions, field solvers, and dno solvers
include("hops_functions.jl")
include("setup_functions.jl")
include("field_tfe_helmholtz_m_and_n.jl")
include("field_tfe_helmholtz_m_and_n_lf.jl")


function two_layer_solve_fast(tau2,zeta_n_m,psi_n_m,
    gamma_u_bar_p,gamma_w_bar_p,N,Nx,f,f_x,pp,alpha_bar,
    gamma_u_bar,gamma_w_bar,Dz,a,b,Nz,M,identy,alpha_bar_p)

  # MSK 7/26/21: Changed the size of U_n_m,W_n_m,ubar_n_m, and wbar_n_m
  # from (Nx,N+1,M+1) to (Nx,M+1,N+1)

  U_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
  W_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
  ubar_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
  wbar_n_m = zeros(Complex{Float64},Nx,M+1,N+1)

  # MSK 7/26/21: Changed the size of G_U_pqrs and J_W_pqrs from
  # (Nx,N+1,N+1,M+1,M+1) to (Nx,M+1,M+1,N+1,N+1)

  G_U_pqrs = zeros(Complex{Float64},Nx,M+1,M+1,N+1,N+1) # G_{p,r}[U_{q,s}]
  J_W_pqrs = zeros(Complex{Float64},Nx,M+1,M+1,N+1,N+1) # J_{p,r}[W_{q,s}]

  #
  # n=0, m=0
  #

  n = 0; m = 0
  Q = zeta_n_m[:,m+1,n+1]
  R = -psi_n_m[:,m+1,n+1]
  U_n_m[:,m+1,n+1],W_n_m[:,m+1,n+1] = AInverse(Q,R,gamma_u_bar_p,gamma_w_bar_p,Nx,tau2)

  # Compute and store G_{p,r}[U_{q,s}]
  q = n; s = m
  # MSK 7/26/21: Changed the size of xi_n_m from (Nx,N-q+1,M-s+1) to (Nx,M-s+1,N-q+1)
  xi_n_m = zeros(Complex{Float64},Nx,M-s+1,N-q+1)
  xi_n_m[:,0+1,0+1] = U_n_m[:,s+1,q+1]
  u_n_m = field_tfe_helmholtz_m_and_n(xi_n_m,f,pp,gamma_u_bar_p,
      alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N-q,M-s,identy,alpha_bar_p)
  G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m,f,pp,Dz,a,Nx,Nz,N-q,M-s)
  for p=0:N-q
    for r=0:M-s
      G_U_pqrs[:,r+1,s+1,p+1,q+1] = G_n_m[:,r+1,p+1]
    end
  end

  # Compute and store J_{p,r}[W_{q,s}]
  q = n; s = m
  # MSK 7/26/21: Changed the size of xi_n_m from (Nx,N-q+1,M-s+1) to (Nx,M-s+1,N-q+1)
  xi_n_m = zeros(Complex{Float64},Nx,M-s+1,N-q+1)
  xi_n_m[:,0+1,0+1] = W_n_m[:,s+1,q+1]
  w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_n_m,f,pp,gamma_w_bar_p,
      alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N-q,M-s,identy,alpha_bar_p)
  J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m,f,pp,Dz,b,Nx,Nz,N-q,M-s)
  for p=0:N-q
    for r=0:M-s
      J_W_pqrs[:,r+1,s+1,p+1,q+1] = J_n_m[:,r+1,p+1]
    end
  end

  #
  # n>0 or m>0
  #
    
  for n=0:N
    for m=0:M
      Q = zeta_n_m[:,m+1,n+1]
      R = -psi_n_m[:,m+1,n+1]
      # MSK 11/01/2021 - Add terms to R from removing the phase
  #     if n > 1 && m > 0
  #       R = R - f_x*1i*alpha_bar.*U_n_m(:,m,n-1) + tau2*f_x*1i*alpha_bar.*W_n_m(:,m,n-1)
  #       if m > 1
  #         R = R - f_x*1i*alpha_bar.*U_n_m(:,m-1,n-1) + tau2*f_x*1i*alpha_bar.*W_n_m(:,m-1,n-1)
  #       end
  #     end
      for r=0:n
        for s=0:m
          if((r<n)||(s<m))
            R = R - G_U_pqrs[:,m-s+1,s+1,n-r+1,r+1] -
                  tau2*J_W_pqrs[:,m-s+1,s+1,n-r+1,r+1]
          end
        end
      end
      if((n>0)||(m>0))
        U_n_m[:,m+1,n+1],W_n_m[:,m+1,n+1] = AInverse(Q,R,gamma_u_bar_p,gamma_w_bar_p,Nx,tau2)
      end
      
      # Compute and store G_{p,r}[U_{q,s}]
      q = n; s = m
      # MSK 7/26/21: Changed the size of xi_n_m from (Nx,N-q+1,M-s+1) to (Nx,M-s+1,N-q+1)
      xi_n_m = zeros(Complex{Float64},Nx,M-s+1,N-q+1)
      xi_n_m[:,0+1,0+1] = U_n_m[:,s+1,q+1]
      u_n_m = field_tfe_helmholtz_m_and_n(xi_n_m,f,pp,gamma_u_bar_p,
          alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N-q,M-s,identy,alpha_bar_p)
      G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m,f,pp,Dz,a,Nx,Nz,N-q,M-s)
      for p=0:N-q
        for r=0:M-s
          G_U_pqrs[:,r+1,s+1,p+1,q+1] = G_n_m[:,r+1,p+1]
        end
      end

      # Compute and store J_{p,r}[W_{q,s}]
      q = n; s = m
      # MSK 7/26/21: Changed the size of xi_n_m from (Nx,N-q+1,M-s+1) to (Nx,M-s+1,N-q+1)
      xi_n_m = zeros(Complex{Float64},Nx,M-s+1,N-q+1)
      xi_n_m[:,0+1,0+1] = W_n_m[:,s+1,q+1]
      w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_n_m,f,pp,gamma_w_bar_p,
          alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N-q,M-s,identy,alpha_bar_p)
      J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m,f,pp,Dz,b,Nx,Nz,N-q,M-s)
      for p=0:N-q
        for r=0:M-s
          J_W_pqrs[:,r+1,s+1,p+1,q+1] = J_n_m[:,r+1,p+1]
        end
      end

    end
  end

  u_n_m = field_tfe_helmholtz_m_and_n(U_n_m,f,pp,gamma_u_bar_p,
      alpha_bar,gamma_u_bar,Dz,a,Nx,Nz,N,M,identy,alpha_bar_p)
  w_n_m = field_tfe_helmholtz_m_and_n_lf(W_n_m,f,pp,gamma_w_bar_p,
      alpha_bar,gamma_w_bar,Dz,b,Nx,Nz,N,M,identy,alpha_bar_p)

  ell_top = 0 + 1
  ell_bottom = Nz + 1

  for n=0:N
    for m=0:M
      ubar_n_m[:,m+1,n+1] = u_n_m[:,ell_top,m+1,n+1]
      wbar_n_m[:,m+1,n+1] = w_n_m[:,ell_bottom,m+1,n+1]
    end
  end

  return U_n_m,W_n_m,ubar_n_m,wbar_n_m
end