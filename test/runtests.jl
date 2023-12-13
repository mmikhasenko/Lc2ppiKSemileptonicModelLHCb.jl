using Test
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecay
using YAML

# get parameters
isobarsinput = YAML.load_file(joinpath(@__DIR__, "..", "data", "particle-definitions.yaml"));
modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "model-definitions.yaml"));
defaultparameters = modelparameters["Default amplitude model"]

# get parameters from json files
# convert to the standard convention
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict=isobarsinput)


# set up the model with a numbering
# 0: Lc, 1:p, 2:pi, 3:K
model = Lc2ppiKModel(; chains, couplings, isobarnames)
@test model isa Lc2ppiKModel

# get a random point in the phase space
σs0 = randomPoint(model.chains[1].tbs.ms)  # (σ1 = m23², σ2 = m31², σ3 = m12²)


# call intensity
_I = unpolarizedintensity(model, σs0)

@test _I isa Real

# call the amplitude
_A = amplitude(model, σs0, [1, 0, 0, 1]) # helici
@assert _A == amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0=1))

@test _A isa Complex
