function f()
    println(global x)
end


function g()
    x = 5
    f()
end


g()