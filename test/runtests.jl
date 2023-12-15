using Test
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecay
using YAML

@testset "Data files are there" begin
    @test isfile(joinpath(@__DIR__, "..", "data", "particle-definitions.yaml"))
    @test isfile(joinpath(@__DIR__, "..", "data", "model-definitions.yaml"))
end

# get parameters
(; particledict, modelparameters) = expose_model_description()

@testset "The model can be are exposed" begin
    @test particledict isa Dict
    @test length(particledict) > 0
    @test modelparameters isa Dict
    @test length(modelparameters) > 0
    @test haskey(modelparameters, "Default amplitude model")
end

defaultparameters = modelparameters["Default amplitude model"]

# get parameters from json files
# convert to the standard convention
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict)

# set up the model with a numbering
# 0: Lc, 1:p, 2:pi, 3:K
model = Lc2ppiKModel(; chains, couplings, isobarnames)


@testset "Pulling default model" begin
    @test model isa Lc2ppiKModel
    @test model == published_model("Default amplitude model")
end

# get a random point in the phase space
σs0 = randomPoint(model.chains[1].tbs.ms)  # (σ1 = m23², σ2 = m31², σ3 = m12²)

# call intensity
_I = unpolarizedintensity(model, σs0)

# call the amplitude
_A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values

@testset "Evaluation of the meeting" begin
    @test _I isa Real
    @test _A isa Complex
    @test _A == amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0=1))
end


@testset "Parameters and couplings" begin
    _names = model.isobarnames
    _HRk = getproperty.(model.chains, :HRk)
    _couplings = model.couplings
    parameter_name(isobarname, recoupling) =
        isobarname * "_{" * string(recoupling.two_λa) * ", " * string(recoupling.two_λb) * "}"

    summary = parameter_name.(_names, _HRk) .=> _couplings
    @test summary isa Lc2ppiKSemileptonicModelLHCb.StaticArrays.SVector{N,T} where {N,T<:Pair{String,ComplexF64}}
end


