using LinearAlgebra
using FFTW
using LinearSolve

"""
    Compute A_nm as a reference for the field solvers.
### Input
    @param alpha_m - ?
    @param f_x - ?
    @param N - ?
    @param M - ?
### Output
    @return the value of A_nm
"""
function A_exp(alpha_m,f_x,N,M)
    Nx = length(f_x);
    A_nm = zeros(Nx,N+1,M+1);
    
    # n = 0, m = 0
    A_nm[:,0+1,0+1] = ones(Nx,1)
    
    # n = 0, m > 0: Do nothing as A_nm[:,0+1,m+1] = 0
    
    # n > 0, m >= 0
    
    for n=1:N
      for m=0:M
        sum = zeros(Nx,1)
        for r=0:m
          sum = sum + A_nm[:,n-1+1,r+1]*(1im*alpha_m[m-r+1])
        end
        A_nm[:,n+1,m+1] = f_x.*sum/n
      end
    end
    
    return A_nm
end

"""
    Compute the inverse of the matrix A.
### Input
    @param Q - ?
### Output
"""
function AInverse(Q,R,gammap,gammapw,Nx,tau2)
    a = zeros(Complex{Float64},Nx,1)
    b = zeros(Complex{Float64},Nx,1)
    Q_hat = fft(Q)
    R_hat = fft(R)

    for j=1:Nx
        det_p = tau2*gammapw[j] + gammap[j]
        a[j] = ((tau2*gammapw[j])*Q_hat[j] + 1im*R_hat[j])/det_p
        b[j] = ((-gammap[j])*Q_hat[j] + 1im*R_hat[j])/det_p
    end

    U = ifft(a)
    W = ifft(b)

    return U,W
end


function chebpts(N)
    return cos.((0:N)*pi/N)
end

"""
Compute and return the Chebyshev differentiation matrix.
Based on Trefethen, Spectral Methods in Matlab, page 54, cheb.m

function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1)."(0:N)'; X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D - diag(sum(D'));
"""    
function cheb(N)
    if N==0
        D=zeros(1,1)
        x=ones(1.)
    else
        x = chebpts(N)
        c = [2; ones(N-1); 2] .* (-1).^(0:N)
        X = repeat(x,outer=[1,N+1])
        dX = X - X'
        D = (c*(1 ./ c)')./(dX + I)
        D = D - diagm(vec(sum(D,dims=2)))
    end
    return D,x
end


function dx(u,p)
    u_x = ifft((1im*p).*fft(u))  
    return u_x
end


function dz(u,Dz,b)
    #u_z = ((2.0/b)*Dz*u')'
    u_z = transpose((2.0/b)*Dz*transpose(u))
    return u_z
end


function E_exp_lf(gamma_m,f,N,M)
    Nx = length(f)
    E_nm_lf = zeros(Complex{Float64},Nx,M+1,N+1)

    # n = 0, m = 0
    E_nm_lf[:,0+1,0+1] = ones(Nx,1)

    # n = 0, m > 0: Do nothing as E_nm_lf(:,0+1,m+1) = 0

    # n > 0, m >= 0

    for n=1:N
        for m=0:M
            sum = zeros(Nx,1)
            for r=0:m
                sum = sum + E_nm_lf[:,r+1,n-1+1]*(1im*gamma_m[m-r+1])
            end
            E_nm_lf[:,m+1,n+1] = -f.*sum/n
        end
    end

    return E_nm_lf
end


function E_exp(gamma_m,f,N,M)
    Nx = length(f)
    E_nm = zeros(Complex{Float64},Nx,M+1,N+1)

    # n = 0, m = 0
    E_nm[:,0+1,0+1] = ones(Nx,1)

    # n = 0, m > 0: Do nothing as E_n_m(:,0+1,m+1) = 0

    # n > 0, m >= 0

    for n=1:N
        for m=0:M
            sum = zeros(Nx,1)
            for r=0:m
                sum = sum + E_nm[:,r+1,n-1+1]*(1im*gamma_m[m-r+1])
            end
            E_nm[:,m+1,n+1] = f.*sum/n
        end
    end

    return E_nm
end


function gamma_exp(alpha_bar,alpha_bar_q,gamma_bar,gamma_bar_q,k_bar,M)
    gamma_q_m = zeros(Complex{Float64},1,M+1)
    gamma_q_m[0+1] = gamma_bar_q

    if(M>=1)
        gamma_q_m[1+1] = (k_bar^2-alpha_bar*alpha_bar_q)/gamma_q_m[0+1]
    end

    if(M>=2)
        gamma_q_m[2+1] = (gamma_bar^2-gamma_q_m[1+1]^2)/(2*gamma_q_m[0+1])
    end

    for m=3:M
        num_sum = 0
        for r=1:m-1
            num_sum = num_sum - gamma_q_m[m-r+1]*gamma_q_m[r+1]
        end
        gamma_q_m[m+1] = num_sum/(2*gamma_q_m[0+1])
    end

    return gamma_q_m
end


function solvebvp_colloc_fast(Uhat,b,alpha,beta,gamma,d_min,n_min,
            r_min,d_max,n_max,r_max,Nx,identy,D,D2,D_start,D_end)

    A = alpha*D2 .+ beta*D .+ reshape(gamma,1,1,Nx).*identy
    A[end,:,:] = repeat(n_min*D_end, outer = [1,1,Nx])
    b[end,:] = r_min
        
    A[1,:,:] = repeat(n_max*D_start, outer = [1,1,Nx])
    A[end,end,:] = A[end,end,:] .+ d_min
    A[1,1,:] = A[1,1,:] .+ d_max
    b[1,:] = r_max

    for j=1:Nx
    #   Apply LinearSolve.jl
    #   utilde = A[:,:,j] \ b[:,j]
    #   Uhat[j,:] = utilde'
        prob = LinearProblem(A[:,:,j],b[:,j])
        utilde = solve(prob)
        Uhat[j,:] = transpose(utilde)
    end

    # TODO - Why doesn't this work?!?
    # b = Matrix{Complex{Float64}}(b)
    # b = reshape(b,(size(b,1),1,size(b,2)))
    # prob = LinearProblem(A,b)
    # soln = solve(prob)

    return Uhat
end


function T_dno(alpha,alphap,gamma,gammap,k2,Nx,M)
    # Tw_dno - Calculates the frequency expansions of gamma_p per the
    # perbutation in delta. See Diffraction Problems: A HOPS/AWE Method
    # section 4.1, (12) - (16).
    #
    # Inputs:
    #
    # alpha - parameter
    # alphap - alpha + tilde(p) = alpha + (2pi / d) * p
    # gamma - parameter
    # gammap - sqrt(k^2 - alpha_p^2) provided p in propogating nodes
    # k2 - alphap(0+1)^2 + gammap(0+1)^2;
    # Nx - number of grid points evaluated on x-axis
    # M - number of iterations for delta pertubation
    #
    # Outputs:
    #
    # T - Result of taylor expansion of gamma_p
    #
    # MSK 4/22/20

    #gammap_m = zeros(Complex{Float64},Nx,M+1)
    # TODO - Will this work?
    gammap_m = zeros(Complex{Float64},Nx,max(M+1,2+1))
    gammap_m[:,0+1] = gammap
    gammap_m[:,1+1] = @. 2*(k2 - alpha*alphap)/(2*gammap)
    gammap_m[:,2+1] = @. (gamma^2 - gammap_m[:,1+1].^2)/(2*gammap)

    for m=3:M
        for r=1:m-1
            gammap_m[:,m+1] = @. gammap_m[:,m+1] -
                (gammap_m[:,m-r+1]*gammap_m[:,r+1])/(2*gammap)
        end
    end

    T = 1im*gammap_m

    return T
end

