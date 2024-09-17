## Task 1
# Returns the sum of even integers between 1 and n inclusive
function sum_of_evens(n::Int)
    sum = 0
    for i in 1:n
        if((i % 2) == 0)
            sum += i
        end
    end

    return sum
end

arg::Int = 5
sum = sum_of_evens(arg)
println("Sum of evens between 1 and $arg (inclusive): $sum" )

## Task 2
# Count the number of occurrences of each integer in an array,
# return as a dictionary
function count_frequencies(arr)
    #arr_size = sizeof(arr) / 4
    freq = Dict{Int, Int}()
    for num in arr
        # Check if the number has already been added to the dictionary
        if haskey(freq, num)
            freq[num] += 1
        else
        # Add the new term to the dictionary
            freq[num] = 1
        end
    end

    return freq
end
arr = [1, 1, 2, 3, 4, 4, 4, 6, 6]
frequencies = count_frequencies(arr)
println(frequencies)

## Task 3
# Overloaded greeting function based on argument type
function greet(name::String)
    return "Hello, $name"
end

function greet(number::Int)
    return "Hello, number $number"
end

function greet(arg)
    return "Hello nerd"
end

println(greet("Doug"))
println(greet(2314))
println(greet(true))


