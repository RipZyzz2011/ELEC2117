using Pkg
Pkg.add("ForwardDiff")
using ForwardDiff

# Evaluate the gradient of a vector function at a point
macro gradient_expr(expr)
    
    quote
        @eval function gradient_function(point::Vector)
            
            #Take the function from the macro expression
            f(x) = $expr
            return ForwardDiff.gradient(f, point)
        end
    end
end
# Define the macro for the expression f(x,y) = x^2 + y^2
# Treat as a vector
@gradient_expr (x[1])^2 + (x[2])^2
point = [1, 2]

grad = gradient_function(point)
println("Gradient at [$(point[1]), $(point[2])]: ",grad)