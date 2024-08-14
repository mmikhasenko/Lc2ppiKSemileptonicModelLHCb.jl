using Lc2ppiKSemileptonicModelLHCb
using Documenter

DocMeta.setdocmeta!(Lc2ppiKSemileptonicModelLHCb, :DocTestSetup, :(using Lc2ppiKSemileptonicModelLHCb); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
  modules = [Lc2ppiKSemileptonicModelLHCb],
  authors = "Mikhail Mikhasenko <mikhail.mikhasenko@cern.ch> and contributors",
  repo = "https://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl/blob/{commit}{path}#{line}",
  sitename = "Lc2ppiKSemileptonicModelLHCb.jl",
  format = Documenter.HTML(; canonical = "https://mmikhasenko.github.io/Lc2ppiKSemileptonicModelLHCb.jl"),
  pages = [
    "index.md"
    [
      file for
      file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
    ]
  ],
)

deploydocs(; repo = "github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl")
