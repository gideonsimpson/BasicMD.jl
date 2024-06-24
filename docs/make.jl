if abspath(PROGRAM_FILE) == @__FILE__
    # When running the `make.jl` file as a script, automatically activate the
    # `docs` environment and dev-install the main package into that environment
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end
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