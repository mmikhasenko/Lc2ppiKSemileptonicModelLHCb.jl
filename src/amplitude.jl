function Lc2ppiKModel(; chains, couplings, isobarnames)

	v = collect(enumerate(isobarnames))
	sort!(v, by = x -> eval(Meta.parse(x[2][3:end-1])))
	sort!(v, by = x -> findfirst(x[2][1], "LDK"))
	order = getindex.(v, 1)

	mypairs = isobarnames[order] .=> zip(couplings[order], chains[order])
	x = Vector{Pair{String, Tuple{Complex, AbstractDecayChain}}}(mypairs)
	# Only p (particle 1) and Lc carry spin; both used reference 1 in
	# the former mixed default (1,2,3,1). The fitted couplings are unchanged.
	ThreeBodyDecay(x; reference_k = 1)
end
