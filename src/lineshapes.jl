
breakup(m², m1², m2²) = sqrtKallenFact(sqrt(m²), sqrt(m1²), sqrt(m2²)) / (2 * sqrt(m²))

function F²(l, p, p0, d)
    pR = p * d
    p0R = p0 * d
    l == 0 && return 1.0
    l == 1 && return (1 + p0R^2) / (1 + pR^2)
    l != 2 && error("l>2 cannot be")
    return (9 + 3p0R^2 + p0R^4) / (9 + 3pR^2 + pR^4)
end

@with_kw struct BuggBreitWigner
    m::Float64
    Γ::Float64
    γ::Float64
end
#
function (BW::BuggBreitWigner)(σ)
    σA = mK^2 - mπ^2 / 2
    @unpack m, Γ, γ = BW
    Γ_dep = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γ_dep)
end


# Flatte1405
@with_kw struct Flatte1405
    m::Float64
    Γ::Float64
    #
    ma::Float64
    mb::Float64
end
#
# Flatte1405(pars::T; kw...) where {T} = Flatte1405(; pars, kw...)
function (BW::Flatte1405)(σ)
    @unpack m, Γ, ma, mb = BW
    p, p0 = breakup(σ, ma^2, mb^2), breakup(m^2, mπ^2, mΣ^2)
    p′, p0′ = breakup(σ, mπ^2, mΣ^2), breakup(m^2, mπ^2, mΣ^2)
    Γ1 = Γ * (p / p0) * m / sqrt(σ)
    Γ2 = Γ * (p′ / p0′) * m / sqrt(σ)
    Γ_tot = Γ1 + Γ2
    1 / (m^2 - σ - 1im * m * Γ_tot)
end

function updatepars(BW::T, pars) where T<:BuggBreitWigner
    fiels = fieldnames(typeof(BW))
    values = [getproperty(BW, f) for f in fiels]
    old_pars = NamedTuple{fiels}(values)
    new_pars = (; old_pars..., pars...)
    return typeof(BW)(; new_pars...)
end

function updatepars(BW::T, pars) where T<:BreitWigner
    fiels = fieldnames(typeof(BW))
    values = [getproperty(BW, f) for f in fiels]
    old_pars = NamedTuple{fiels}(values)
    new_pars = (; old_pars..., pars...)
    return typeof(BW)(; new_pars...)
end
