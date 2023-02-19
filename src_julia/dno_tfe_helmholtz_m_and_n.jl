using FFTW

function dno_tfe_helmholtz_m_and_n(unm,f,p,Dz,a,Nx,Nz,N,M)

    Gnm = zeros(Complex{Float64},Nx,M+1,N+1)
    
    ell_bottom = Nz + 1
    f_x = ifft( (1im*p).*fft(f) )
    
    for n=0:N
      for m=0:M
        u_z = dz(unm[:,:,m+1,n+1],Dz,a)
        Gnm[:,m+1,n+1] = -u_z[:,ell_bottom]
        if(n>=1)
          u_x = dx(unm[:,:,m+1,n-1+1],p)
          Gnm[:,m+1,n+1] = Gnm[:,m+1,n+1] + f_x.*u_x[:,ell_bottom]
        
          Gnm[:,m+1,n+1] = Gnm[:,m+1,n+1] + (1.0/a)*(f.*Gnm[:,m+1,n-1+1])
        end
        if(n>=2)
          u_x = dx(unm[:,:,m+1,n-2+1],p)
          Gnm[:,m+1,n+1] = Gnm[:,m+1,n+1] - (1.0/a)*(f.*(f_x.*u_x[:,ell_bottom]))
    
          u_z = dz(unm[:,:,m+1,n-2+1],Dz,a)
          Gnm[:,m+1,n+1] = Gnm[:,m+1,n+1] - f_x.*(f_x.*u_z[:,ell_bottom])
        end
      end
    end
    
    return Gnm
  end