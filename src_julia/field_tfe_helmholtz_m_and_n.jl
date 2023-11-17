using FFTW
using LinearAlgebra

"""
Compute the approximate solution (unm) in the upper field.
### Input
- `xi_lf_n_m` -- a tensor representing the partial solution in the upper field
- `f` -- a test function representing the grating surface
- `p` -- an integer where tilde_p = (2π/d)p and d is the periodicity of the grating interface
- `gammap` -- a numerical constant in the upper field for all wave numbers p
- `alpha` -- a numerical constant
- `gamma` -- a numerical constant in the upper field 
- `Dz` -- the partial derivative with respect to the z component
- `a` -- the artificial boundary imposed at the top of the upper layer
- `Nx` -- the number of discretization points
- `Nz` -- the number of collocation points
- `N` -- the maximum number of Taylor orders for the interfacial perturbation
- `M` -- the maximum number of Taylor orders for the frequency perturbation
- `identy` -- the identity matrix
- `alphap` --  a numerical constant at all wave numbers p
### Output
- `unm` --  a tensor representing the approximate solution in the upper field
"""
function field_tfe_helmholtz_m_and_n(xi_n_m,f,p,gammap,
                alpha,gamma,Dz,a,Nx,Nz,N,M,identy,alphap)

  unm = zeros(Complex{Float64},Nx,Nz+1,M+1,N+1)

  #k2 = p[0+1]^2 + gammap[0+1]^2
  k2 = alphap[0+1]^2 + gammap[0+1]^2

  ell_top = 0 + 1
  xi_n_m_hat = fft(xi_n_m,(1,))
  f_x = real(ifft( (1im*p).*fft(f) ))

  ll = (0:Nz)
  z_min = 0; z_max = a
  # MSK 7/22/21 - Parameters for solvebvp_colloc
  D = (2/(z_max-z_min))*Dz
  D2 = D*D
  D_start = D[1,:]
  D_end = D[end,:]
  # End MSK 7/22/21
  tilde_z = @. cos(pi*ll/Nz)
  z = @. ((z_max-z_min)/2.0)*(tilde_z - 1.0) + z_max

  f_full = repeat(f, outer=[1,Nz+1])
  f_x_full = repeat(f_x, outer=[1,Nz+1])
  # TODO - Speedup the line below?
  a_minus_z_full = repeat(a .- transpose(z), outer=[Nx,1])

  Uhat = zeros(Complex{Float64},Nx,Nz+1)

  Tu = T_dno(alpha,p,gamma,gammap,k2,Nx,M)

  # n=0 and m=0

  for ell=0:Nz
    unm[:,ell+1,0+1,0+1] = ifft( exp.(1im*gammap*z[ell+1]).*xi_n_m_hat[:,0+1,0+1] )
  end

  A1_xx = -(2.0/a)*f_full
  A1_xz = @. -(1.0/a)*(a_minus_z_full)*f_x_full
  A1_zx = A1_xz
  #A1_zz = 0
    
  A2_xx = (1.0/a^2)*f_full.^2
  A2_xz = @. (1.0/a^2)*(a_minus_z_full)*(f_full*f_x_full)
  A2_zx = A2_xz
  A2_zz = @. (1.0/a^2)*((a_minus_z_full)^2)*(f_x_full^2)
    
  B1_x = (1.0/a)*f_x_full
  #B1_z = 0
    
  B2_x = @. -(1.0/a^2)*f_full*f_x_full
  B2_z = @. -(1.0/a^2)*(a_minus_z_full)*(f_x_full^2)

  S1 = -(2.0/a)*f_full
  S2 = @. (1.0/a^2)*f_full^2

  for n=0:N
    for m=0:M
      
      # Form Fnm, Jnm
      Fnm = zeros(Complex{Float64},Nx,Nz+1)
      Jnm = zeros(Complex{Float64},Nx,1)
      
      if(n>=1)
        u_x = dx(unm[:,:,m+1,n-1+1],p)
        temp = @. A1_xx*u_x
        Fnm = Fnm - dx(temp,p)
        temp = @. A1_zx*u_x
        Fnm = Fnm - dz(temp,Dz,a)
        temp = @. B1_x*u_x
        Fnm = Fnm - temp
      
        u_z = dz(unm[:,:,m+1,n-1+1],Dz,a)
        temp = @. A1_xz*u_z
        Fnm = Fnm - dx(temp,p)
        #A1_zz = 0
        #B1_z = 0
      
        temp = @. 2*1im*alpha*S1*u_x
        Fnm = Fnm - temp
        temp = @. (gamma^2)*S1*unm[:,:,m+1,n-1+1]
        Fnm = Fnm - temp
      end
    
      if(m>=1)
        u_x = dx(unm[:,:,m-1+1,n+1],p)
        temp = @. 2*1im*alpha*u_x
        Fnm  = Fnm - temp
        temp = @. 2*(gamma^2)*unm[:,:,m-1+1,n+1]
        Fnm  = Fnm - temp
      end
    
      if(n>=1 && m>=1)
        u_x = dx(unm[:,:,m-1+1,n-1+1],p)
        temp = @. 2*1im*alpha*S1*u_x
        Fnm  = Fnm - temp
        temp = @. 2*(gamma^2)*S1*unm[:,:,m-1+1,n-1+1]
        Fnm  = Fnm - temp
      end
    
      if(n>=2)
        u_x = dx(unm[:,:,m+1,n-2+1],p)
        temp = @. A2_xx*u_x
        Fnm = Fnm - dx(temp,p)
        temp = @. A2_zx*u_x
        Fnm = Fnm - dz(temp,Dz,a)
        temp = @. B2_x*u_x
        Fnm = Fnm - temp
        
        u_z = dz(unm[:,:,m+1,n-2+1],Dz,a)
        temp = @. A2_xz*u_z
        Fnm = Fnm - dx(temp,p)
        temp = @. A2_zz*u_z
        Fnm = Fnm - dz(temp,Dz,a)
        temp = @. B2_z*u_z
        Fnm = Fnm - temp
        
        temp = @. 2*1im*alpha*S2*u_x
        Fnm = Fnm - temp
        temp = @. (gamma^2)*S2*unm[:,:,m+1,n-2+1]
        Fnm = Fnm - temp
      end
      
      if(m>=2)
        temp = @. (gamma^2)*unm[:,:,m-2+1,n+1]
        Fnm = Fnm - temp
      end
    
      if(n>=1 && m>=2)
        temp = @. (gamma^2)*S1*unm[:,:,m-2+1,n-1+1]
        Fnm = Fnm - temp
      end
    
      if(n>=2 && m>=1)
        u_x = dx(unm[:,:,m-1+1,n-2+1],p)
        temp = @. 2*1im*alpha*S2*u_x
        Fnm = Fnm - temp
        temp = @. 2*(gamma^2)*S2*unm[:,:,m-1+1,n-2+1]
        Fnm = Fnm - temp
      end
    
      if(n>=2 && m>=2)
        temp = @. (gamma^2)*S2*unm[:,:,m-2+1,n-2+1]
        Fnm = Fnm - temp
      end
      
      for r=0:m-1
        Jnm = Jnm + ifft( (Tu[:,m-r+1]).*fft(unm[:,ell_top,r+1,n+1], (1,)), (1,))
      end
      if(n>=1)
        for r=0:m
          Snm = ifft( (Tu[:,m-r+1]).*fft(unm[:,ell_top,r+1,n-1+1], (1,)), (1,) )
          Jnm = Jnm - (1.0/a)*f.*Snm
        end
      end
      
      # Solve elliptic equation
      
      Fnmhat = fft(Fnm,(1,))
      Jnmhat = fft(Jnm,(1,))
      
      b = transpose(Fnmhat)
      alphaalpha = 1.0
      betabeta = 0.0
      #gammagamma = gamma*gamma - p.^2 - 2*alpha.*p
      #gammagamma = gammap.^2
      gammagamma = @. k2 - alphap^2
      d_min = 1.0
      n_min = 0.0
      r_min = xi_n_m_hat[:,m+1,n+1]
      d_max = -1im*gammap
      n_max = 1.0
      r_max = Jnmhat
      
      Uhat = solvebvp_colloc_fast(Uhat,b,alphaalpha,betabeta,gammagamma,
          d_min,n_min,r_min,d_max,n_max,r_max,Nx,identy,D,D2,D_start,D_end)
      
      if((n>0)||(m>0))
        unm[:,:,m+1,n+1]=ifft(Uhat, (1,))
      end
      
    end
  end

  return unm

end