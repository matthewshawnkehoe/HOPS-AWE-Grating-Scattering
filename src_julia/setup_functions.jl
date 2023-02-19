# setup_functions.jl
# Four functions used to setup information in other programs

function setup_2d(Nx,L,alpha,beta)
    xx = (L/Nx)*(0:Nx-1)
    kk = (2.0*pi/L)*[0:Nx/2-1;-Nx/2:-1]
      
    alphap = alpha .+ kk
    betap = zeros(Complex{Float64},Nx,1)
    
    kappa = sqrt(alpha^2 + beta^2)
    
    value = kappa^2 .- alphap.^2
    rho = abs.(value)
    theta = angle.(value)
    
    # 'angle' returns -pi < theta < pi
    # If theta < 0 add 2*pi so that 0 < theta < 2*pi
    
    for j=1:Nx
      if(theta[j]<0)
        theta[j] = theta[j] + 2*pi
      end
    end
    
    # Always choose betap so that Im{betap} > 0
    # Thus exp(i*betap*y) bounded for y -> infty
    # Thus exp(-i*betap*y) bounded for y -> -infty
    
    for j=1:Nx
      betap[j] = sqrt(rho[j])*exp(1im*theta[j]/2.0)
    end

    # TODO - Necessary?
    # for j=1:Nx
    #   val = sqrt(rho[j])*exp(1im*theta[j]/2.0)
    #   if real(val) < 1e-15
    #     betap[j] = 1im * imag(val)
    #   else
    #     betap[j] = val
    #   end
    # end
    
    eep = exp.(1im*alpha*xx)
    eem = exp.(-1im*alpha*xx)

    return xx,kk,alphap,betap,eep,eem
end


function setup_xi_u_nu_u_n_m(A,r,xx,pp,alpha_bar_p,
                gamma_bar_p,f,f_x,Nx,N,M)
    alpha_bar = alpha_bar_p[0+1]
    gamma_bar = gamma_bar_p[0+1]
    k_bar = sqrt(alpha_bar^2 + gamma_bar^2)
    
    pp_r = pp[r+1]
    alpha_bar_r = alpha_bar_p[r+1]
    gamma_bar_r = gamma_bar_p[r+1]
    
    # AWE approximation of gamma_q
    
    gamma_r_m = gamma_exp(alpha_bar,alpha_bar_r,
             gamma_bar,gamma_bar_r,k_bar,M)
    
    # HOPS/AWE approximation of xi_q and nu_q
    
    E_n_m = E_exp(gamma_r_m,f,N,M)
    upper =  A*exp.(1im*pp_r*xx)
    
    xi_r_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    nu_r_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    
    for n=0:N
      for m=0:M
        # TODO - Why does E_n_m[:,m+1,n+1] need to be transposed?
        xi_r_n_m[:,m+1,n+1] = upper.*E_n_m[:,m+1,n+1]
        for ell=0:m
          nu_r_n_m[:,m+1,n+1] = nu_r_n_m[:,m+1,n+1] +
              (-1im*gamma_r_m[m-ell+1]).*xi_r_n_m[:,ell+1,n+1]
        end
        if(n>0)
          nu_r_n_m[:,m+1,n+1] = nu_r_n_m[:,m+1,n+1] +
              f_x.*(1im*pp_r).*xi_r_n_m[:,m+1,n-1+1]
        end
      end
    end
    
    return xi_r_n_m,nu_r_n_m
end


function setup_xi_w_nu_w_n_m(A,r,xx,pp, alpha_bar_p,
            gamma_bar_p,f,f_x,Nx,N,M)
    alpha_bar = alpha_bar_p[0+1]
    gamma_bar = gamma_bar_p[0+1]
    k_bar = sqrt(alpha_bar^2 + gamma_bar^2)
    
    pp_r = pp[r+1]
    alpha_bar_r = alpha_bar_p[r+1]
    gamma_bar_r = gamma_bar_p[r+1]
    
    # AWE approximation of gamma_q
    
    gamma_r_m = gamma_exp(alpha_bar,alpha_bar_r,
        gamma_bar,gamma_bar_r,k_bar,M)
    
    # HOPS/AWE approximation of xi_q and nu_q
    
    E_n_m = E_exp_lf(gamma_r_m,f,N,M)
    lower = A*exp.(1im*pp_r*xx)
    
    xi_r_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    nu_r_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    
    for n=0:N
      for m=0:M
        xi_r_n_m[:,m+1,n+1] = lower.*E_n_m[:,m+1,n+1]
        for ell=0:m
          nu_r_n_m[:,m+1,n+1] = nu_r_n_m[:,m+1,n+1] -
              lower.*(1im*gamma_r_m[m-ell+1]).*E_n_m[:,ell+1,n+1]
        end
        if(n>0)
          nu_r_n_m[:,m+1,n+1] = nu_r_n_m[:,m+1,n+1] -
              lower.*f_x.*(1im*pp_r).*E_n_m[:,m+1,n-1+1]
        end    
      end
    end
    
    return xi_r_n_m,nu_r_n_m
end


function setup_zeta_psi_n_m(xx,pp,alpha_bar,gamma_u_bar,
                f,f_x,Nx,N,M)
    gamma_u_m = zeros(1,M+1)
    gamma_u_m[0+1] = gamma_u_bar
    gamma_u_m[1+1] = gamma_u_bar
    
    alpha_u_m = zeros(1,M+1)
    alpha_u_m[0+1] = alpha_bar
    alpha_u_m[1+1] = alpha_bar
    
    E_n_m = E_exp_lf(gamma_u_m,f,N,M)
    
    zeta_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    psi_n_m = zeros(Complex{Float64},Nx,M+1,N+1)
    
    # MSK 02/05/23 Initialize ell to zero to account for the case where
    # m = 0 and n > 0
    ell = 0
    
    for n=0:N
      for m=0:M
        zeta_n_m[:,m+1,n+1] = -E_n_m[:,m+1,n+1]
        for ell=0:m
          psi_n_m[:,m+1,n+1] = psi_n_m[:,m+1,n+1] +
              (1im*gamma_u_m[m-ell+1]).*E_n_m[:,ell+1,n+1]
        end

        if(n>0)
          psi_n_m[:,m+1,n+1] = psi_n_m[:,m+1,n+1] +
              f_x.*(1im*alpha_u_m[m-ell+1]).*E_n_m[:,ell+1,n-1+1]
        end    
      end
    end
    
    return zeta_n_m,psi_n_m
end