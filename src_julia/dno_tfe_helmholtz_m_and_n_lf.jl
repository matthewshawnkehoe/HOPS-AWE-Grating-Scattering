using FFTW

"""
Compute the lower layer DNO.
### Input
- `wmn` -- a tensor representing the solution in the lower layer
- `f` -- a test function representing the grating surface
- `p` -- an integer where tilde_p = (2π/d)p and d is the periodicity of the grating interface
- `Dz` -- the partial derivative with respect to the z component
- `b` -- the artificial boundary imposed at the bottom of the lower layer
- `Nx` -- the number of discretization points
- `Nz` -- the number of collocation points
- `N` --- the maximum number of Taylor orders for the interfacial perturbation
- `M` -- the maximum number of Taylor orders for the frequency perturbation
### Output
- `Jnm` -- a tensor representing the lower layer DNO
"""
function dno_tfe_helmholtz_m_and_n_lf(wnm,f,p,Dz,b,Nx,Nz,N,M)

  Jnm = zeros(Complex{Float64},Nx,M+1,N+1)

  ell_top = 0 + 1
  f_x = ifft( (1im*p).*fft(f) )

  for n=0:N
    for m=0:M
      w_z = dz(wnm[:,:,m+1,n+1],Dz,b)
      Jnm[:,m+1,n+1] = w_z[:,ell_top]
      if(n>=1)
        w_x = dx(wnm[:,:,m+1,n-1+1],p)
        Jnm[:,m+1,n+1] = Jnm[:,m+1,n+1] - f_x.*w_x[:,ell_top]
      
        Jnm[:,m+1,n+1] = Jnm[:,m+1,n+1] - (1.0/b)*(f.*Jnm[:,m+1,n-1+1])
      end
      if(n>=2)
        w_x = dx(wnm[:,:,m+1,n-2+1],p)
        Jnm[:,m+1,n+1] = Jnm[:,m+1,n+1] - (1.0/b)*(f.*(f_x.*w_x[:,ell_top]))

        w_z = dz(wnm[:,:,m+1,n-2+1],Dz,b)
        Jnm[:,m+1,n+1] = Jnm[:,m+1,n+1] + f_x.*(f_x.*w_z[:,ell_top])
      end
    end
  end

  return Jnm
end