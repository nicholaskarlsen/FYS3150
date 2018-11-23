using BenchmarkTools

N = 1e6
x = zeros(convert(Int, N))

function f(array)
    for val in array
        val = Threads.threadid()
    end
end

@btime f(x)