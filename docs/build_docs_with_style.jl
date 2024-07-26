CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

 import Pkg
 Pkg.pkg"add Documenter@1"

 # Load documenter
 using Documenter
 using DocumenterTools: Themes
 ENV["JULIA_DEBUG"] = "Documenter"

 # For easier debugging when downloading from a specific branch.
 github_user = "JuliaDynamics"
 branch = "master"
 download_path = "https://raw.githubusercontent.com/$github_user/doctheme/$branch"

 import Downloads
 for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
     Downloads.download("$download_path/$file", joinpath(@__DIR__, file))
 end

 # create the themes
 for w in ("light", "dark")
     header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
     theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
     write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
 end

 # compile the themes
 Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
 Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

 # Download and apply CairoMakie plotting style
 using CairoMakie
 Downloads.download("$download_path/style.jl", joinpath(@__DIR__, "style.jl"))
 include("style.jl")

 using DocumenterInterLinks

links = InterLinks(
    "ComplexityMeasures" => (
        "https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/",
        "https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/objects.inv",
        joinpath(@__DIR__, "inventories", "Documenter.toml")
    ),
    "StateSpaceSets" => (
        "https://juliadynamics.github.io/DynamicalSystemsDocs.jl/statespacesets/stable/",
        "https://juliadynamics.github.io/DynamicalSystemsDocs.jl/statespacesets/stable/objects.inv",
        joinpath(@__DIR__, "inventories", "Documenter.toml")
    ),
);

 function build_docs_with_style(pages, modules...; bib = nothing, authors = "George Datseris", draft = false, kwargs...)
     settings = (
         modules = [modules...],
         format = Documenter.HTML(
             prettyurls = CI,
             assets = [
                 asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
             ],
             collapselevel = 3,
             size_threshold = 1250 * 2^10, 
             size_threshold_warn = 750 * 2^10,
             example_size_threshold = 100 * 2^10,
             size_threshold_ignore = [
                "examples/examples_associations.md", 
                "examples/examples_independence.md",
                "examples/examples_infer_graphs.md",
             ]
         ),
         sitename = "$(modules[1]).jl",
         authors,
         pages,
         draft,
         doctest = false,
         warnonly = false,
         checkdocs = :exported,
         kwargs...
     )

     if isnothing(bib)
         makedocs(; settings...)
     else
         makedocs(; plugins=[bib, links], settings...)
     end

     if CI
         deploydocs(
             repo = "github.com/JuliaDynamics/$(modules[1]).jl.git",
             target = "build",
             push_preview = true
         )
     end

 end 