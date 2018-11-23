using BenchmarkTools

x = range(0, stop=convert(Int64, 1e6), length=convert(Int64, 1e6))

function f(array)
    for val in array
        val+=1
    end
end

function g(array)
    Threads.@threads for val in array
        val+=1
    end
end



println("Number of threads: ", Threads.nthreads())
print("Without paralellization")
@btime f(x)
print("With paralellization")
@btime g(x)