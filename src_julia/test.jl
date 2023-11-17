# Padding - Method 1
# lpad(string, n, "p"): Make a string at least n columns wide when printed, by padding on the left with copies of p.
display(lpad(4,7,"0"))
# similar for rpad
display(rpad(4,3,"0"))

# Padding - Method 2
using Images

# https://juliaimages.org/stable/function_reference/#ImageFiltering.padarray

x = [1 2 ;3 4]
display(x)
y = padarray(x, Pad(1, 1))

# Padding - Method 3
# https://github.com/JuliaArrays/PaddedViews.jl

# PaddedViews provides a simple wrapper type, PaddedView, to add "virtual" padding to any array without copying data. 
# Edge values not specified by the array are assigned a fillvalue. Multiple arrays may be "promoted" to have common indices using the paddedviews function.

# PaddedView arrays are read-only, meaning that you cannot assign values to them. The original array may be extracted using A = parent(P), 
# where P is a PaddedView.

a = collect(reshape(1:9, 3, 3))
display(PaddedView(-1, a, (3,6)))
PaddedView(-1, a, (1:5,1:5), (2:4,2:4))
# A = parent(a)
# display(A)

# range or LinRange works for a comparable version of np.linspace()