#=
These calls were required to make julia work on my computer, similar may be required
if you can not run main.py. See the README file
=#

ENV["PYTHON"] = "/home/nicholas/anaconda2/bin/python"
using Pkg
Pkg.build("PyCall")