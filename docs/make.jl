push!(LOAD_PATH,"../src/")
using BasicMD
using Documenter
makedocs(
         sitename = "BasicMD.jl",
         modules  = [BasicMD],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/gideonsimpson/BasicMD.jl",
)