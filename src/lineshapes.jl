
breakup(m², m1², m2²) = sqrtKallenFact(sqrt(m²), sqrt(m1²), sqrt(m2²)) / (2 * sqrt(m²))

@with_kw struct BuggBreitWigner
    m::Float64
    Γ::Float64
    γ::Float64
end
#
function (BW::BuggBreitWigner)(σ)
    mK = 0.493677
    mπ = 0.13957018
    σA = mK^2 - mπ^2 / 2
    @unpack m, Γ, γ = BW
    Γ_dep = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    1 / (m^2 - σ - 1im * m * Γ_dep)
end

@with_kw struct BuggBreitWignerExpFF
    m::Float64
    Γ::Float64
    γ::Float64
    α::Float64
end
#
function (BW::BuggBreitWignerExpFF)(σ)
    mK = 0.493677
    mπ = 0.13957018
    σA = mK^2 - mπ^2 / 2
    @unpack m, Γ, γ, α = BW
    Γ_dep = (σ - σA) / (m^2 - σA) * Γ * exp(-γ * σ)
    q = breakup(σ, mK^2, mπ^2)
    exp(-α * q^2) / (m^2 - σ - 1im * m * Γ_dep)
end

# Flatte1405
@with_kw struct Flatte1405
    m::Float64
    Γ1::Float64
    Γ2::Float64
    #
    ma::Float64
    mb::Float64
end
#
function (BW::Flatte1405)(σ)
    mπ = 0.13957018
    mΣ = 1.18937
    @unpack m, Γ1, Γ2, ma, mb = BW
    # 
    p_norm = breakup(m^2, mπ^2, mΣ^2)
    p = breakup(σ, ma^2, mb^2) 
    p′ = breakup(σ, mπ^2, mΣ^2)
    # 
    gsq1 = Γ1 / (2p_norm) * m
    gsq2 = Γ2 / (2p_norm) * m
    Γ1_dep = gsq1 * 2p / sqrt(σ) 
    Γ2_dep = gsq2 * 2p′/ sqrt(σ)
    Γ_dep = Γ1_dep + Γ2_dep
    1 / (m^2 - σ - 1im * m * Γ_dep)
end

# updatepars

function updatepars(BW::T, pars) where T<:Union{BreitWigner, BuggBreitWigner, BuggBreitWignerExpFF, Flatte1405}
    fiels = fieldnames(typeof(BW))
    values = [getproperty(BW, f) for f in fiels]
    old_pars = NamedTuple{fiels}(values)
    new_pars = (; old_pars..., pars...)
    return typeof(BW)(; new_pars...)
end
