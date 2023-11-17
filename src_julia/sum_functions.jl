# MSK 01/29/2023 - Added conditional statements for padesum/taylorsum when the coeff matrix is of size one

using ToeplitzMatrices
using Polynomials
using DSP

function fcn_sum(SumType,f_n_m,Eps,delta,Nx,N,M)

    f = zeros(Complex{Float64},Nx,1)
    coeff = zeros(Complex{Float64},N+1,M+1)

    for j=1:Nx
      for r=0:N
        for s=0:M
          coeff[r+1,s+1] = f_n_m[j,s+1,r+1]
        end
      end
      if(SumType==1)
        # Change to taylorsum_2_coeff if M,N >= 6
        f[j] = taylorsum_2_coeff(coeff,Eps,delta,N,M)[1]
      elseif(SumType==2)
        f[j] = padesum2(coeff,Eps,delta,N,M)[1]
      elseif(SumType==3)
        f[j] = padesum2_safe(coeff,Eps,delta,N,M)[1]
      end
    end

    return f
end

function fcn_sum_fast(SumType,f_n_m,Eps,delta,Nx,N,M)
    f = zeros(Complex{Float64},Nx,1)
    
    for j=1:Nx
      if(SumType==1)
        # Change to taylorsum_2_coeff if M,N >= 6
        f[j] = taylorsum_2_coeff(f_n_m[:,:,j],Eps,delta,N,M)[1]
      elseif(SumType==2)
        f[j] = padesum2(f_n_m[:,:,j],Eps,delta,N,M)[1]
      else 
        f[j] = padesum2_safe(f_n_m[:,:,j],Eps,delta,N,M)[1]
      end
    end
    
    return f
end


function padesum_safe(c,Eps,M)
    # padesum - Sums a truncated Taylor series via Pade approximation.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_0,...,c_N]
    # Eps - Value at which to sum
    # M - Numerator and denominator degree (M = N/2)
    #
    # Outputs:
    #
    # psum - Pade approximant evaluated at Eps
    # a - Numerator coefficients: [a_0,...,a_M]
    # b - Denominator coefficients: [1,b_1,...,b_M]
    #
    # Note: Assumes that N = 2*M (otherwise ignores c_N)
    #
    # DPN 2/7/12
    # DPN 11/25/14 (safe version)

    N_true = -1
    for j=1:2*M+1
      if(abs(c[j])>1e-14)
        N_true = j
      end
    end
    M_safe = floor(Int, N_true/2.0)

    if(M_safe==0)
      a = c[1]
      b = [1]
    else
      # MSK 01/29/2023 replaced Matlab toeplitz with Julia Toeplitz
      H = Toeplitz(c[M_safe+1:2*M_safe-1+1],c[M_safe+1:-1:1+1])
      bb = -H\c[M_safe+1+1:2*M_safe+1]
      b = [1;bb]
      aa = conv(b,c[0+1:M_safe+1])
      a = aa[0+1:M_safe+1]
    end

    if length(c)==1
      psum = polyval(a[1],Eps)/polyval(b[1],Eps)
    else
      psum = polyval(a,Eps)./polyval(b,Eps)
    end

    return psum,a,b
end


function padesum(c,Eps,M)
    # padesum - Sums a truncated Taylor series via Pade approximation.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_0,...,c_N]
    # Eps - Value at which to sum
    # M - Numerator and denominator degree (M = N/2)
    #
    # Outputs:
    #
    # psum - Pade approximant evaluated at Eps
    # a - Numerator coefficients: [a_0,...,a_M]
    # b - Denominator coefficients: [1,b_1,...,b_M]
    #
    # Note: Assumes that N = 2*M (otherwise ignores c_N)
    #
    # DPN 2/7/12
      
    if(M==0)
      a = c[1]
      b = [1]
    else
      # MSK 01/29/2023 replaced Matlab toeplitz with Julia Toeplitz
      H = Toeplitz(c[M+1:2*M-1+1],c[M+1:-1:1+1])
      bb = -H\c[M+1+1:2*M+1]
      b = [1;bb]
      aa = conv(b,c[0+1:M+1])
      a = aa[0+1:M+1]
    end

    if length(c)==1
      psum = polyval(a[1],Eps)/polyval(b[1],Eps)
    else
      psum = polyval(a,Eps)./polyval(b,Eps)
    end

    return psum,a,b
end


function padesum2_safe(c,Eps,delta,N,M)
    # padesum2 - Uses Pade approximation to sum a truncated Taylor series.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_{0,0},...,c_{0,N},...,c_{M,0},...,c_{M,N}]
    # Eps, delta - Value at which to sum
    # N,M - Degree of truncated Taylor series
    #
    # Outputs:
    #
    # psum - Taylor sum evaluated at Eps, delta
    #
    # DPN 12/20/15
    # MSK 7/21/21 Optimized double for loop
      
    rho = sqrt(Eps^2 + delta^2)
    # TODO - Should this be angle ?
    theta = atan(delta,Eps)
      
    # Form the ctilde
    N_M_min = min(N,M)
      
    coeff = zeros(Complex{Float64},N_M_min+1,1)
    c1 = cos(theta).^(0:N_M_min)
    c2 = sin(theta).^(0:N_M_min)
    for p=0:N_M_min
      for q=0:p
        #ii = index_nm(p-q,q,N,M);
        coeff[p+1] = coeff[p+1] + c[p-q+1,q+1]*c1[p-q+1]*c2[q+1]
      end
    end
      
    N_M_min_over_2 = floor(Int, N_M_min/2.0)

    if length(c) == 1
      psum = padesum_safe(coeff[1],rho,N_M_min_over_2)[1]
    else
      psum = padesum_safe(coeff,rho,N_M_min_over_2)
    end
    
    return psum[1]
end


function padesum2(c,Eps,delta,N,M)
    # padesum2 - Uses Pade approximation to sum a truncated Taylor series.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_{0,0},...,c_{0,N},...,c_{M,0},...,c_{M,N}]
    # Eps, delta - Value at which to sum
    # N,M - Degree of truncated Taylor series
    #
    # Outputs:
    #
    # psum - Taylor sum evaluated at Eps, delta
    #
    # DPN 12/20/15
    # MSK 7/21/21 Optimized double for loop
      
    rho = sqrt(Eps^2 + delta^2)
    # TODO - Should this be angle ?
    theta = atan(delta,Eps)
      
    # Form the ctilde
      
    N_M_min = min(N,M)
      
    coeff = zeros(Complex{Float64},N_M_min+1,1)
    c1 = @. cos(theta)^(0:N_M_min)
    c2 = @. sin(theta)^(0:N_M_min)
    for p=0:N_M_min
      for q=0:p
        coeff[p+1] = coeff[p+1] + c[p-q+1,q+1]*c1[p-q+1]*c2[q+1]
      end
    end
      
    N_M_min_over_2 = floor(Int, N_M_min/2.0)

    if length(c) == 1
      psum = padesum(coeff[1],rho,N_M_min_over_2)[1]
    else
      psum = padesum(coeff,rho,N_M_min_over_2)
    end
    
    return psum[1]
end


function taylorsum_2_coeff(c,Eps,delta,N,M)
    # taylorsum - Sums a truncated Taylor series.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_0,...,c_N]
    # Eps - Value at which to sum for N
    # Delta - Value at which to sum for M
    # N - Degree of truncated Taylor series
    # M - Degree of truncated Taylor series
    #
    # Outputs:
    #
    # tsum - Taylor sum evaluated at Eps and Delta
    #
    # DPN 2/7/12
    # MSK 7/21/21 Optimized double for loop

    rho = sqrt(Eps^2 + delta^2)
    # TODO - Should this be angle ?
    theta = atan(delta,Eps)

    # Form the ctilde

    N_M_min = min(N,M)

    coeff = zeros(Complex{Float64},N_M_min+1,1)
    c1 = @. cos(theta)^(0:N_M_min)
    c2 = @. sin(theta)^(0:N_M_min)
    for p=0:N_M_min
      for q=0:p
        #ii = index_nm(p-q,q,N,M)
        coeff[p+1] = coeff[p+1] + c[p-q+1,q+1]*c1[p-q+1]*c2[q+1]
      end
    end
    
    N_M_min_over_2 = floor(Int, N_M_min/2.0)

    if length(c) == 1
      tsum = taylorsum(coeff[1],rho,N_M_min_over_2)[1]
    else
      tsum = taylorsum(coeff,rho,N_M_min_over_2)
    end

    return tsum
end


function taylorsum(c,Eps,N)
    # taylorsum - Sums a truncated Taylor series.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_0,...,c_N]
    # Eps - Value at which to sum
    # N - Degree of truncated Taylor series (unused!)
    #
    # Outputs:
    #
    # tsum - Taylor sum evaluated at Eps
    #
    # DPN 2/7/12

    #tsum = polyval(wrev(c),Eps)
    #tsum = polyval(wrev(c(1:N+1)),Eps)
    if length(c) == 1
      tsum = polyval(c[1],Eps)
    else
      tsum = polyval(c,Eps)
    end

    #p = Poly(c)
    # tsum = polyval(c, Eps)

    if length(c) == 1
      return tsum[1]
    else
      return tsum
    end
end


""" Evaluate a polynomial vs. Horner's rule """
function polyval(c,x)
    # based on Matlab's Horner's rule 
    m = length(x)
    y = zeros(Complex{Float64},m)
    
    if length(c) > 0
        y[:] .= c[end]
    end
    for i=2:length(c)
        y = @. x*y + c[end-i+1]
    end
    return y
end


function taylorsum2(c,Eps,delta,N,M)
    # taylorsum2 - Sums a truncated Taylor series.
    #
    # Inputs:
    #
    # c - Taylor series coefficients: [c_{0,0},...,c_{0,N},...,c_{M,0},...,c_{M,N}]
    # Eps, delta - Value at which to sum
    # N,M - Degree of truncated Taylor series
    #
    # Outputs:
    #
    # tsum - Taylor sum evaluated at Eps, delta
    #
    # MSK 7/11/2021

    # TODO - Rewrite differently?
    vander = @. (Eps^(0:N))*(delta^(0:M)')
    # TODO - Should this be sum(c.*vander)?
    tsum = sum(sum(c.*vander))

    # Alternative Implementation (slightly faster for large M,N)
    # c1 = zeros(N+1,1)
    # c2 = zeros(1,M+1)
    # c1[1,1]=1
    # c2[1,1]=1
    # c1[2:end] = cumprod(repmat(Eps,N,1))
    # c2[2:end] = cumprod(repmat(delta,1,M))
    # tsum = sum(c.*c1.*c2, 'all')

    if length(c) == 1
      return tsum[1]
    else
      return tsum
    end
end


function vol_fcn_sum(SumType,u_n_m,Eps,delta,Nx,Nz,N,M)
    u = zeros(Complex{Float64},Nx,Nz+1)
    coeff = zeros(Complex{Float64},N+1,M+1)

    for j=1:Nx
      for ell=0:Nz
        for r=0:N
          for s=0:M
            coeff[r+1,s+1] = u_n_m[j,ell+1,s+1,r+1]
          end
        end
        if(SumType==1)
          u[j,ell+1] = taylorsum2(coeff,Eps,delta,N,M)[1]
        elseif(SumType==2)
          u[j,ell+1] = padesum2(coeff,Eps,delta,N,M)[1]
        elseif(SumType==3)
          u[j,ell+1] = padesum2_safe(coeff,Eps,delta,N,M)[1]
        end
      end
    end

    return u
end