using Test
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
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

@testset "Exact values of amplitude and intensity" begin
    @test _I ≈ 6302.587073290807
    @test _A ≈ -0.7072336123636882 - 14.887651207373576im
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


@testset "MeasuredParameter with cases" begin
    @test MeasuredParameter(1.0, 0.1, 0.2, 0.3) isa MeasuredParameter
    @test MeasuredParameter(1.0, 0.1) isa MeasuredParameter
    @test MeasuredParameter("0.91 ± 0.04 ± 0.35 ± 0.04") ==
          MeasuredParameter(0.91, 0.04, 0.35, 0.04)
    @test MeasuredParameter("0.91 ± 0.04") ==
          MeasuredParameter(0.91, 0.04, 0.0, 0.0)
end


@testset "Alternative models" begin
    list_of_models = keys(modelparameters) |> collect
    # exclude the last one, which is the LS model
    # it cannot be loaded, not implemented
    list_of_models_but_few = filter(list_of_models) do x
        x != "Alternative amplitude model obtained using LS couplings" &&
            x != "Alternative amplitude model with an additional overall exponential form factor exp(-alpha q^2) multiplying Bugg lineshapes. The exponential parameter is indicated as ``alpha''" &&
            x != "Alternative amplitude model with free radial parameter d for the Lc resonance, indicated as dLc" &&
            x != "Alternative amplitude model with free L(1405) Flatt'e widths, indicated as G1 (pK channel) and G2 (Sigmapi)"
    end
    @test map(list_of_models_but_few) do modelname
        published_model(modelname) isa Lc2ppiKModel
    end |> all
end
