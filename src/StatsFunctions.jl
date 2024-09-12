# Returns the mean value of an array
function mean_value(data)
    sum = 0
    for i in data
        sum += i
    end

    return sum / length(data)
end
# Returns the variance of an array
function variance_value(data)
    mean = mean_value(data)
    var = 0
    for i in data
        var += (i - mean)^2
    end

    var = var / (length(data))
    return var
end

#Compute the median of an array
function median_value(data)
    # Obtain the length of the array
    sorted = sort(data)
    len = length(data)
    median = 0
    # Check if it is even or odd
    if((len % 2) != 0)
        median = sorted[div((len + 1), 2)]
    else
        median = (sorted[(div(len,2) + 1)] + sorted[div(len,2)]) / 2
    end

    return median
end

function standard_deviation(data)
    return sqrt(variance_value(data))
end

function z_score(data, x)
    sigma = standard_deviation(data)
    mu = mean_value(data)
    return (x - mu) / sigma
end

        

