using CommunityEcology
using Documenter

DocMeta.setdocmeta!(CommunityEcology, :DocTestSetup, :(using CommunityEcology); recursive=true)

makedocs(;
    modules=[CommunityEcology],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/CommunityEcology.jl/blob/{commit}{path}#{line}",
    sitename="CommunityEcology.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/CommunityEcology.jl",
        repolink="https://github.com/EvoArt/CommunityEcology.jl",
    ),
    pages=["Home" => "index.md"],
)

deploydocs(repo="github.com/EvoArt/CommunityEcology.jl")
