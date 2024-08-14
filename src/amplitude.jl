function Lc2ppiKModel(; chains, couplings, isobarnames)

    v = collect(enumerate(isobarnames))
    sort!(v, by = x -> eval(Meta.parse(x[2][3:end-1])))
    sort!(v, by = x -> findfirst(x[2][1], "LDK"))
    order = getindex.(v, 1)

    mypairs = isobarnames[order] .=> zip(couplings[order], chains[order])
    x = Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(mypairs)
    ThreeBodyDecay(x)
end
