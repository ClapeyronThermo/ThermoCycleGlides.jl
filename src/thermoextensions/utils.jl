

function isentropic_compressor(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray{TZ}, fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, TZ<:Real}
    TT = promote_type(T1, T2, T3, T4, TZ)
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in, z,phase = :vapour)::TT
    # T_in = Tproperty(fluid,p_in,h_in,z,enthalpy,phase = :vapour)
    T_isen_out = Tproperty(fluid,p_out,s_isen,z,entropy,phase = :vapour)::TT
    h_isen = enthalpy(fluid,p_out,T_isen_out,z,phase = :vapour)::TT

    # h_isen = Clapeyron.PS.enthalpy(fluid, p_out, s_isen, z)

    ha =  (h_in + ((h_isen - h_in)/η_isen))::TT
    #  T_out = Clapeyron.PH.temperature(fluid,p_out,ha,z)
    Tcrit,pcrit,_ = crit_mix(fluid,z)
    if p_out < pcrit
        T_dew = dew_temperature(fluid,p_out,z)[1]::TT
        h_dew = enthalpy(fluid,p_out,T_dew,z,phase=:vapour)::TT
        if ha < h_dew
        # @warn "Fixing outlet of compressor at saturation temperature"
            return enthalpy(fluid,p_out,T_dew,z,phase = :vapour)::TT
        end
    end

    return ha
end

function isentropic_pump(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray, fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in, z,phase = :liquid)
    T_isen_out = Tproperty(fluid,p_out,s_isen,z,entropy,phase = :liquid)
    h_isen = enthalpy(fluid,p_out,T_isen_out,z,phase = :liquid)
    return h_in + (h_isen - h_in) / η_isen
end

function isentropic_expander(p_in::T1,p_out::T2,η_isen::T3,h_in::T4,z::AbstractVector,fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in,z) 
    T_isen_out = Tproperty(fluid,p_out,s_isen,z,entropy,phase = :vapour)
    h_isen = enthalpy(fluid,p_out,T_isen_out,z,phase = :vapour)
    h_out = h_in - (h_in - h_isen) * η_isen
    # force outlet to be gaseous 
    dt = Clapeyron.dew_temperature(fluid, p_out, z)[1]
    h_gas = Clapeyron.enthalpy(fluid, p_out, dt, z,phase = :gas)
    # @show   values(dt), h_gas, h_out, p_in, p_out
    if h_out < h_gas
        # @warn "The outlet enthalpy is below the gas enthalpy. Adjusting to gas phase."
        h_out = h_gas 
    end

    return h_out
end

function IHEX_Q(fluid::EoSModel,ϵ::T1,T_in_left::T2,p_in_left::T3,T_in_right::T4,p_in_right::T5,z::AbstractVector{T}) where {T<:Real,T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
    T_in_left <= T_in_right && @warn  "Function assumse left fluid is hotter. This didn't happen. Now heat transfer will be negative."
    c_in_left = Clapeyron.isobaric_heat_capacity(fluid,p_in_left,T_in_left,z)
    c_in_right = Clapeyron.isobaric_heat_capacity(fluid,p_in_right,T_in_right,z)
    C = min(c_in_left,c_in_right)
    Qmax = C*(T_in_left - T_in_right)
    Q = Qmax*ϵ
    return Q
end


function glide_match_coeff(prob::ThermoCycleGlides.ThermoCycleProblem,sol::SolutionState;N::Int = 20)

end


function compute_isentropic_exponent(fluid::EoSModel,p,T,z)
    Cp = isobaric_heat_capacity(fluid,p,T,z)
    Cv = isochoric_heat_capacity(fluid,p,T,z)
    return Cp/Cv
end

"""
https://data.dtu.dk/articles/dataset/Numerical_models_for_the_design_and_analysis_of_heat_pumps_with_zeotropic_mixtures/6825443?file=13709117
"""
function off_design_compressor_relation(fluid::EoSModel,z,η_isen_design,built_in_volume_ratio;p_ref = 101325.0, T_ref = 300)
    @assert 0.0 <= η_isen_design <= 1 "design point isentropic efficiecny should be between (0,1)"
    @assert dew_temperature(fluid,p_ref,z)[1] ≤ T_ref "Reference state should be gas. Change p_ref or T_ref" 
    κ = compute_isentropic_exponent(fluid,p_ref,T_ref,z)
    πi = built_in_volume_ratio^κ
    γ = (κ - 1)/κ
    f(x) = η_isen_design*((x)^(γ) - 1)/(πi^γ - (γ*πi^(-1/κ))*(πi - x) - 1)
    return f
end

export off_design_compressor_relation, compute_isentropic_exponent
