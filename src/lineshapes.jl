
breakup(m², m1², m2²) = sqrtKallenFact(sqrt(m²), sqrt(m1²), sqrt(m2²)) / (2 * sqrt(m²))

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
    # 
    p_norm = breakup(m^2, mπ^2, mΣ^2)
    p = breakup(σ, ma^2, mb^2) 
    p′ = breakup(σ, mπ^2, mΣ^2)
    # 
    gsq = Γ / (2p_norm) * m
    Γ1 = gsq * 2p / sqrt(σ) 
    Γ2 = gsq * 2p′/ sqrt(σ)
    Γ_tot = Γ1 + Γ2
    1 / (m^2 - σ - 1im * m * Γ_tot)
end

# updatepars

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
