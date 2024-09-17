using Pkg
Pkg.add("Zytgote")
using Zytgote

# Take a single variable expressiona and evaluate its value and derivative at a point
macro diff_macro(expr)
    quote
        @eval function value_func(x)
            f(x) = $expr
            return f(x)
        end
        @eval function deriv_func(x)
            f(x) = $expr
            deriv = Zytgote.derivative(f, x)
            return deriv
        end
    end
end

@diff_macro x^2 + 3*x + 5
point = 2
value = value_func(point)
derivative = deriv_func(point)

println("Value of x^2 + 3x + 5 at x = $point: $value")
println("Derivative of x^2 + 3x + 5 at x = $point: $derivative")
