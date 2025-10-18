using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using Lc2ppiKSemileptonicModelLHCb
import Lc2ppiKSemileptonicModelLHCb: lineshape_mismatch
using YAML
using Test

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
    @test model isa ThreeBodyDecay
    @test model == published_model("Default amplitude model")
end

# get a random point in the phase space
σs0 = Invariants(model.chains[1].tbs.ms;
    σ1=0.7980703453578917,
    σ2=3.6486261122281745)

# call intensity
_I = unpolarized_intensity(model, σs0)

# call the amplitude
_A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values


@testset "Evaluation of the meeting" begin
    @test _I isa Real
    @test _A isa Complex
    @test _A == amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0=1))
end



@testset "Exact values of amplitude and intensity" begin
    @test _I ≈ 9345.853380852352
    @test _A ≈ -45.1323269502508 + 54.85942516648639im
    # 
    # @test model.chains[1].Xlineshape(σs0.σ2) ≈
    #       model.chains[2].Xlineshape(σs0.σ2) ≈ -0.5636481410171861 + 0.13763637759224928im
    # # 
    # @test model.chains[21].Xlineshape(σs0.σ1) ≈
    #       model.chains[22].Xlineshape(σs0.σ1) ≈
    #       model.chains[23].Xlineshape(σs0.σ1) ≈
    #       model.chains[24].Xlineshape(σs0.σ1) ≈ 2.1687201455088894 + 23.58225917009096im

end

@testset "Parameters and couplings" begin
    _names = model.names
    _HRk = getproperty.(model.chains, :HRk)
    _couplings = model.couplings
    parameter_name(isobarname, vertex) =
        isobarname * "_{" * string(vertex.h.two_λa) * ", " * string(vertex.h.two_λb) * "}"

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
# let
    list_of_models = keys(modelparameters) |> collect
    # exclude the last one, which is the LS model
    # it cannot be loaded, not implemented
    list_of_models_but_few = filter(list_of_models) do x
        x != "Alternative amplitude model obtained using LS couplings" &&
            # requires new lineshape
            x != "Alternative amplitude model with an additional overall exponential form factor exp(-alpha q^2) multiplying Bugg lineshapes. The exponential parameter is indicated as ``alpha''" &&
            # requires reading d and replacing
            x != "Alternative amplitude model with free radial parameter d for the Lc resonance, indicated as dLc" &&
            # requires reading flatte coupkings
            x != "Alternative amplitude model with free L(1405) Flatt'e widths, indicated as G1 (pK channel) and G2 (Sigmapi)"
    end
    @test map(list_of_models_but_few) do modelname
        # @show modelname
        published_model(modelname) isa ThreeBodyDecay
    end |> all
end
