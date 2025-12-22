# Lc2pKpi Semileptonic Model from LHCb publication

[![10.1007/JHEP07(2023)228](<https://zenodo.org/badge/doi/10.1007/JHEP07(2023)228.svg>)](<https://doi.org/10.1007/JHEP07(2023)228>)
[![GPLv3+ license](https://img.shields.io/badge/License-GPLv3+-blue.svg)](https://www.gnu.org/licenses/gpl-3.0-standalone.html)

## Overview

`Lc2ppiKSemileptonicModelLHCb.jl` is a Julia package developed for processing and analyzing $\Lambda^+_c \to p \pi^+ K^-$ decays with LHCb experiment in the semileptonic production of Lc. This package is separates Julia code from the [`polarimetry`](https://github.com/ComPWA/polarimetry) project, enabling to set up and run the model quickly.
Helicity couplings and other parameter values are taken from a recent study by the LHCb Collaboration[^1] and its [supplementary material](https://cds.cern.ch/record/2824328/files).

[^1]: Amplitude analysis of the $\Lambda^+_c \to p K^- \pi^+$ decay and $\Lambda^+_c$ baryon polarization measurement in semileptonic beauty hadron decays (2022) [[link]](https://inspirehep.net/literature/2132745)

## Installation

To install `Lc2ppiKSemileptonicModelLHCb.jl`, use the Julia package manager:

```julia
using Pkg
Pkg.add("https://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl")  # this code
```

### New project and environment

In case you install a fresh julia on your laptop,
follow the steps to create a new project folder, install `Lc2ppiKSemileptonicModelLHCb` to dependences.

1. start julia in a terminal
2. locate yourself with `pwd()`
3. go to a new project folder with `cd("path")` and `mkdir("folder")` if needed
4. once you are in a project folder, do

```julia
] activate
```

5. check that you have a clean environment with

```julia
] st
```

6. add dependences

```julia
]
add http://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl
```

## Usage

After installation, you can import the package and begin your analysis:

```julia
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays

model = published_model("Default amplitude model")

# module is a simple combination of `couplings` and `chains` arrays
# where the chain is rather flat structure of decay information
model.chains[3] |> dump

# get a random point in the phase space
σs0 = randomPoint(masses(model))  # (σ1 = m23², σ2 = m31², σ3 = m12²)
two_λs0 = ThreeBodySpins(1,0,0; two_h0=1) # helicity values

# call intensity
_I = unpolarized_intensity(model, σs0)

# call the amplitude
_A = amplitude(model, σs0, two_λs0)  # pars: model, mandelstam variables, helicity values

# take TBS algebra for dalitz plot
σs_test = let m2Kπ = 0.79, m2πp = 3.64
    Invatriants(masses(model); σ1=m2Kπ, σ2=m2πp)
end
#
# evaluate what you want
unpolarized_intensity(model, σs_test) # full model
amplitude(model, σs0, two_λs0)
amplitude(model.chains[2], σs0, two_λs0)  # for just 1 chain, number 2
```

## Contributing

Contributions to `Lc2ppiKSemileptonicModelLHCb.jl` are welcome. Just create an issue.

## Reference to LHCb Analysis

This project is based on the original analysis conducted by the LHCb collaboration. For detailed information, refer to the [LHCb experiment publications](https://lhcb-public.web.cern.ch/en/Publications/en).

## Related Projects

- **Polarimetry**: Visit the [`polarimetry`](https://github.com/ComPWA/polarimetry) repository for more information on the comprehensive framework for polarization analysis.
- **Dalitz Plot Decomposition**: the model contention are aligned with [Dalitz-plot decomposition](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.034033).
  See [ThreeBodyDecays.jl](https://github.com/mmikhasenko/ThreeBodyDecays.jl) and [SymbolicThreeBodyDecays.jl](https://github.com/mmikhasenko/SymbolicThreeBodyDecays.jl) for further details.

## License

`Lc2ppiKSemileptonicModelLHCb.jl` is available under [GPLv3+ license](https://github.com/mmikhasenko/Lc2ppiKSemileptonicModelLHCb.jl/blob/main/LICENSE).
