
function readjson(path)
    f = read(path, String)
    return JSON.parse(f)
end

function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end

function definechaininputs(key, dict; tbs, parities)
    @unpack mass, lineshape, l = dict
    #
    k = Dict('K' => 1, 'D' => 3, 'L' => 2)[first(key)]
    #
    jp_R = str2jp(dict["jp"])
    parity = jp_R.p
    two_j = jp_R.two_j
    #
    i, j = ij_from_k(k)
    #
    @unpack two_js = tbs
    #
    reaction_ij = jp_R => (SpinParity(two_js[i], parities[i]), SpinParity(two_js[j], parities[j]))
    reaction_Rk(P0) = SpinParity(two_js[4], P0) => (jp_R, SpinParity(two_js[k], parities[k]))
    #
    two_LS = vcat(possible_ls.(reaction_Rk.(('+', '-')))...)
    minLS = first(sort(vcat(two_LS...); by=x -> x[1]))
    minL = div(minLS[1], 2)
    #
    Hij = VertexFunction(ParityRecoupling(two_js[i], two_js[j], reaction_ij), BlattWeisskopf{l}(1.5))
    # 
    Xlineshape = build_lineshape(lineshape |> Symbol |> eval, dict)
    return (; k, Xlineshape, Hij, two_j, parity, minL)
end

function build_lineshape(::Type{BreitWigner}, dict)
    m = dict["mass"] / 1000
    Γ = dict["width"] / 1000
    ma = dict["ma"] / 1000
    mb = dict["mb"] / 1000
    l = dict["l"]
    d = dict["d"]
    return BreitWigner(; m, Γ, ma, mb, l, d)
end
function build_lineshape(::Type{Flatte1405}, dict)
    m = dict["mass"] / 1000
    Γ1 = dict["G1"] / 1000
    Γ2 = dict["G2"] / 1000
    ma = dict["ma"] / 1000
    mb = dict["mb"] / 1000
    return Flatte1405(; m, Γ1, Γ2, ma, mb)
end
function build_lineshape(::Type{BuggBreitWigner}, dict)
    m = dict["mass"] / 1000
    Γ = dict["width"] / 1000
    γ = dict["gamma"]
    return BuggBreitWigner(; m, Γ, γ)
end
function build_lineshape(::Type{BuggBreitWignerExpFF}, dict)
    m = dict["mass"] / 1000
    Γ = dict["width"] / 1000
    γ = dict["gamma"]
    α = dict["alpha"]
    return BuggBreitWignerExpFF(; m, Γ, γ, α)
end
# shape parameters
function parseshapedparameter(par_name)
    keytemp = r"([M,G]|G[12]|gamma|alpha)"
    nametemp = r"([L,K,D]\([0-9]*\))"
    m = match(keytemp * nametemp, par_name)
    m === nothing && error("The name of the shared parameter, $(par_name), is not recognized!")
    return (key=m[1], isobarname=m[2])
end

function keyname2symbol(key)
    key == "M" && return :m
    key == "G" && return :Γ
    key == "G1" && return :Γ1
    key == "G2" && return :Γ2
    key == "gamma" && return :γ
    key == "alpha" && return :α
    error("The name of the shared parameter, $(key), is not recognized!")
end

function replacementpair(parname, val)
    @unpack key, isobarname = parseshapedparameter(parname)
    s = keyname2symbol(key)
    v = MeasuredParameter(val).val
    isobarname => eval(:(NamedTuple{($(QuoteNode(s)),)}($(v))))
end


function build_reaction_kinematics(reaction, particledict)
    initial_state = reaction["initial-state"]
    final_state = reaction["final-state"]

    # Extract masses (convert MeV to GeV)
    m0 = particledict[initial_state]["mass"] / 1000
    m1 = particledict[final_state[1]]["mass"] / 1000
    m2 = particledict[final_state[2]]["mass"] / 1000
    m3 = particledict[final_state[3]]["mass"] / 1000

    ms = ThreeBodyMasses(m1=m1, m2=m2, m3=m3, m0=m0)

    # Extract spins and build parities array
    jp_0 = str2jp(particledict[initial_state]["jp"])
    jp_1 = str2jp(particledict[final_state[1]]["jp"])
    jp_2 = str2jp(particledict[final_state[2]]["jp"])
    jp_3 = str2jp(particledict[final_state[3]]["jp"])

    parities = [jp_1.p, jp_2.p, jp_3.p, '±']  # Force '±' for parity-violating decay

    tbs = ThreeBodySystem(ms, ThreeBodySpins(jp_1.two_j, jp_2.two_j, jp_3.two_j; two_h0=jp_0.two_j))

    return (; tbs, parities)
end

function parse_model_dictionaries(modeldict; particledict)
    # Extract reaction information and build kinematics
    reaction = modeldict["reaction"]
    (; tbs, parities) = build_reaction_kinematics(reaction, particledict)

    # 1) get isobars - NOW PASS tbs and parities
    isobars = Dict()
    for (key, lineshape) in modeldict["lineshapes"]
        dict = Dict{String,Any}(particledict[key])
        dict["lineshape"] = lineshape
        isobars[key] = definechaininputs(key, dict; tbs, parities)  # Pass tbs and parities
    end

    # 3) get parameters
    defaultparameters = modeldict["parameters"]
    shapeparameters = filter(x -> x[1] != 'A', keys(defaultparameters))
    parameterupdates = [
        replacementpair(parname, defaultparameters[parname])
        for parname in shapeparameters]

    for (p, u) in parameterupdates
        BW = isobars[p].Xlineshape
        isobars[p] = merge(isobars[p],
            (Xlineshape=updatepars(BW, u),))
    end

    # 3) get couplings
    couplingkeys = collect(filter(x -> x[1:2] == "Ar", keys(defaultparameters)))
    isobarnames = map(x -> x[3:end-1], couplingkeys)

    terms = []
    for parname in couplingkeys
        c_re_key = "Ar" * parname[3:end]
        c_im_key = "Ai" * parname[3:end]
        value_re = MeasuredParameter(defaultparameters[c_re_key]).val
        value_im = MeasuredParameter(defaultparameters[c_im_key]).val
        value = value_re + 1im * value_im
        #
        c0, d = parname2decaychain(parname, isobars; tbs)
        #
        push!(terms, (c0 * value, d))
    end

    chains = getindex.(terms, 2)
    couplings = getindex.(terms, 1)

    couplings .*= lineshape_mismatch.(chains)

    (; chains, couplings, isobarnames)
end



"""
    expose_model_description()

Reads the data, returns disctionaries of isobars and list of a dictionary of models with their parameters.
"""
function expose_model_description()
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "model-definitions.yaml"))
    return (; modelparameters, particledict)
end


"""
    published_model(modelname="Default amplitude model")

Reads the data, performs coupling conversion, returns the model with the parameters.
"""
function published_model(modelname="Default amplitude model")
    # 
    (; modelparameters, particledict) = expose_model_description()
    defaultparameters = modelparameters[modelname]
    (; chains, couplings, isobarnames) =
        parse_model_dictionaries(defaultparameters; particledict)
    model = Lc2ppiKModel(; chains, couplings, isobarnames)
    # 
    return model
end