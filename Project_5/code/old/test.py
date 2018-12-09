from julia import Main

Main.include("test.jl")
x = Main.eval("hello()")

print x