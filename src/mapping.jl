
function selectindexmap(isobarname)
    #
    couplingindexmap = Dict(
        r"[L|D].*" => Dict(
            '1' => (1, 0),
            '2' => (-1, 0)),
        r"K\(892|1410\)" => Dict(
            '1' => (0, 1),
            '2' => (-2, 1),
            '3' => (2, -1),
            '4' => (0, -1)),
        r"K\(700|1430\)" => Dict(
            '1' => (0, -1),
            '2' => (0, 1)))
    #
    m = filter(keys(couplingindexmap)) do k
        match(k, isobarname) !== nothing
    end
    length(m) == 0 && error("no match found for $(isobarname)!")
    return couplingindexmap[first(m)]
end

function parname2decaychain(parname, isobars)
    isobarname = parname[3:end-1]
    #
    @unpack k, Hij, two_j, Xlineshape, parity, minL = isobars[isobarname]
    #
    couplingindex = parname[end]
    two_λR, two_λk = selectindexmap(isobarname)[couplingindex]
    two_λR′, two_λk′, c′ =
        couplingLHCb2DPD(two_λR, two_λk; k, two_j, parity)
    HRk = VertexFunction(NoRecoupling(two_λR′, two_λk′), BlattWeisskopf{minL}(5.0))
    (c′, DecayChain(; k, Xlineshape, Hij, HRk, two_j, tbs))
end

function couplingLHCb2DPD(two_λR, two_λk; k, parity, two_j)
    if k == 2
        @assert two_λk == 0
        c′ = -(2 * (parity == '+') - 1) / sqrt(two_j + 1)
        return (-two_λR, two_λk, c′)
    elseif k == 3
        c′ = -(2 * (parity == '+') - 1) * minusone()^(two_j // 2 - 1 // 2) / sqrt(two_j + 1)
        return (-two_λR, two_λk, c′)
    end
    k != 1 && error("cannot be!")
    c′ = 1.0 / sqrt(two_j + 1)
    return (two_λR, -two_λk, c′)
end

lineshape_mismatch(dc::DecayChain{<:Flatte1405}) = 1.0
lineshape_mismatch(dc::DecayChain{<:BuggBreitWignerMinL}) = 1.0

function lineshape_mismatch(dc::DecayChain)
    minL = orbital_momentum(dc.HRk.ff)
    # 
    dR, dΛc = 1.5, 5.0 # /GeV
    @unpack pars, l = dc.Xlineshape
    @unpack m = pars

    ms = masses(dc)
    @unpack k = dc
    i,j = ij_from_k(k)
    mi, mj, mk = ms[i], ms[j], ms[k]
    m0 = ms[4]
    # 
    p0 = breakup(m^2, mi^2, mj^2)
    q0 = breakup(m0^2, m^2, mk^2)
    # 
    factor_l = (l == 2) ? 9 : 1
    factor_minL = (minL == 2) ? 9 : 1
    1 / (dR * p0)^l * sqrt(F²(l, 0, p0, dR) * factor_l) * 
        1 / (dΛc * q0)^minL * sqrt(F²(minL, 0, q0, dΛc) * factor_minL)
end

"""
The relation is
```math
A^{DPD}_{λ₀,λ₁} = (-1)^{½-λ₁} A^{LHCb}_{λ₀,-λ₁}
```
"""
amplitudeLHCb2DPD(A) =
    [A[1, 2] -A[1, 1]
        A[2, 2] -A[2, 1]]
#
