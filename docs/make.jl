push!(LOAD_PATH,"../src/")
using BasicMD
using Documenter
makedocs(checkdocs=:none,
         sitename = "BasicMD.jl",
         modules  = [BasicMD],
         pages=[
                "Home" => "index.md",
                "Sampling"=>"sample1.md",
                "Samplers"=>["samplers/metropolis1.md", 
                "samplers/nonmetropolis1.md"],
                "Utilities"=>["utils/opts1.md"],
                "Examples" => ["examples/sample_traj1.md", 
                "examples/sample_obs1.md", "examples/sample_con1.md"]
               ])
deploydocs(;
    repo="github.com/gideonsimpson/BasicMD.jl",
)