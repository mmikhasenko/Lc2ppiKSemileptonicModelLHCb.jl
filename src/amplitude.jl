struct Lc2ppiKModel{N, D <: DecayChain, L <: Number}
	chains::SVector{N, D}
	couplings::SVector{N, L}
	isobarnames::SVector{N, <:String}
end

function Lc2ppiKModel(; chains, couplings, isobarnames)
	N = length(chains)
	N != length(couplings) && error("Length of couplings does not match the length of the chains")
	N != length(isobarnames) && error("Length of isobarnames does not match the length of the chains")
	#
	v = collect(zip(chains, couplings, isobarnames))
	sort!(v, by = x -> eval(Meta.parse(x[3][3:end-1])))
	sort!(v, by = x -> findfirst(x[3][1], "LDK"))
	#
	sort_chains, sort_couplings, sort_isobarnames =
		getindex.(v, 1), getindex.(v, 2), getindex.(v, 3)
	#
	sv_sort_chains = (SVector{N, <:DecayChain})(sort_chains)
	sv_sort_couplings = SVector{N}(sort_couplings)
	sv_sort_isobarnames = SVector{N}(sort_isobarnames)
	#
	Lc2ppiKModel(sv_sort_chains, sv_sort_couplings, sv_sort_isobarnames)
end


amplitude(model::Lc2ppiKModel, σs, two_λs) =
	sum(c * amplitude(d, σs, two_λs) for (c, d) in zip(model.couplings, model.chains))

unpolarizedintensity(model::Lc2ppiKModel, σs) =
	sum(abs2, amplitude(model, σs, two_λs)
			  for two_λs in itr(tbs.two_js))
# 
masses(model::Lc2ppiKModel) = first(model.chains).tbs.ms
