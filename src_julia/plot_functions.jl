using Plots; pythonplot()
# using Plots
# gr()
# using PyCall
# pygui(:tk)
# using PyPlot
default(show = true)
using LaTeXStrings
using Measures


function fourier_repr_lipschitz(P,x)
    f = 0
    for k = 1:ceil(P/2)
        f = @. f + 8/(pi^2*(2*k-1)^2) * cos((2*k-1)*x)
    end
    
    f_x = 0
    for k = 1:ceil(P/2)
        f_x = @. f_x + (-8*(2*k-1))/(pi^2*(2*k-1)^2) * sin((2*k-1)*x)
    end

    g1 = plot(x,f,label=L"$f(x)$")
    xaxis!(L"$x$")
    yaxis!(L"$f(x)$")
    title!(L"Lipschitz profile for $f$")
    g2 = plot(x,f_x,label=L"$f_x(x)$")
    xaxis!(L"$x$")
    yaxis!(L"$f_x(x)$")
    title!(L"Lipschitz profile for $f_x$")
    fig = plot(g1,g2,layout=(1,2),margin=2mm,dpi=200)
    display(fig)
    
    return f,f_x
end

function fourier_repr_rough(P,x)
    f = 0;
    for k = 1:P
        f = @. f + 96*(2*k^2*pi^2 - 21)/(125*k^8) * cos(k*x)
    end

    f_x = 0
    for k = 1:P
        f_x = @. f_x + (-96*k*(2*k^2*pi^2 - 21))/(125*k^8) * sin(k*x)
    end

    g1 = plot(x,f,label=L"$f(x)$")
    xaxis!(L"$x$")
    yaxis!(L"$f(x)$")
    title!(L"Lipschitz profile for $f$")
    g2 = plot(x,f_x,label=L"$f_x(x)$")
    xaxis!(L"$x$")
    yaxis!(L"$f_x(x)$")
    title!(L"Lipschitz profile for $f_x$")
    fig = plot(g1,g2,layout=(1,2),margin=2mm,dpi=200)
    display(fig)

    return f,f_x
end


# Make error plots
function plot_errors(sum_type,N,M,err_G,err_U,err_ubar,err_J,err_W,err_wbar)

    p1 = contourf([0:N],[0:M],log10.(err_G)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in G ($sum_type)")
    p2 = contourf([0:N],[0:M],log10.(err_U)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in U ($sum_type)")
    p3 = contourf([0:N],[0:M],log10.(err_ubar)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in U_bar ($sum_type)")
    p4 = contourf([0:N],[0:M],log10.(err_J)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in J ($sum_type)")
    p5 = contourf([0:N],[0:M],log10.(err_W)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in W ($sum_type)")
    p6 = contourf([0:N],[0:M],log10.(err_wbar)')
    xaxis!("n")
    yaxis!("m")
    title!("Error in W_bar ($sum_type)")

    return p1,p2,p3,p4,p5,p6
end


# Make MMS error plots
function plot_mms_error(omega, Eps, err_pade)
    c = contourf(omega,Eps,log10.(abs.(err_pade)),color=:hot, reuse=false)
    xaxis!(L"$\omega$")
    yaxis!(L"$\varepsilon$")
    title!("Relative Error")
    return c
end


# Plot the energy defect (D) for Taylor sumation
function plot_energy_defect_taylor(c_ed,PlotLambda,omega,lambda,Eps,ee_taylor)
    if(PlotLambda==0)
        return contourf!(c_ed,omega,Eps,log10.(abs.(ee_taylor)))
    else
        return contourf!(c_ed,lambda,Eps,log10.(abs.(ee_taylor)))
    end
end


# Plot the energy defect (D) for Pade sumation
function plot_energy_defect_pade(c_ed,PlotLambda,omega,lambda,Eps,ee_pade)
    if(PlotLambda==0)
        return contourf!(c_ed,omega,Eps,log10.(abs.(ee_pade)))
    else
        return contourf!(c_ed,lambda,Eps,log10.(abs.(ee_pade)))
    end
end


# Plot the reflectivity map (R) for Taylor summation
function plot_refl_map_taylor(c_rm,PlotRelative,omega,lambda,Eps,ru_taylor,ru_flat)
    if(PlotRelative==0)
        RR = ru_taylor
    else
        RR = ru_taylor./ru_flat
    end
    if(PlotLambda==0)
        contourf!(c_rm,omega,Eps,RR)
    else
        contourf!(c_rm,lambda,Eps,RR)
    end
end


# Plot the reflectivity map (R) for Pade summation
function plot_refl_map_pade(c_rm,PlotRelative,omega,lambda,Eps,ru_pade,ru_flat)
    if(PlotRelative==0)
        RR = ru_pade
    else
        RR = ru_pade./ru_flat
    end
    if(PlotLambda==0)
        contourf!(c_rm,omega,Eps,RR)
    else
        contourf!(c_rm,lambda,Eps,RR)
    end
end