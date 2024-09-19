using Pkg
Pkg.add("ForwardDiff")
using ForwardDiff

# Evaluate the gradient of a vector function at a point
macro gradient_expr(expr)
    
    return quote
        @eval function gradient_function(point::Vector)
            
            #Take the function from the macro expression
            f = $(expr)
            return ForwardDiff.gradient(x -> f(x...), point)
        end
    end
end
# Define the macro for the expression f(x,y) = x^2 + y^2
# Treat as a vector
@gradient_expr (x,y) -> x^2 + y^2
point = [1, 2]
gradient_function = gradient_function
grad = gradient_function(point)
println("Gradient at [$(point[1]), $(point[2])]: ",grad)