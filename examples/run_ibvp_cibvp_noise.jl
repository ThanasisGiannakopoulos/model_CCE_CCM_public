
using model_CCE_CCM
using Parameters
using Random

@with_kw struct noise_ibvp <: IBVP
    noise_amp_ϕ1  :: Float64 = 0.0
    noise_amp_ψv1 :: Float64 = 0.0
    noise_amp_ψ1  :: Float64 = 0.0
    # speeds (absolute value)
    vϕ1 :: Float64 = 1.0 # right moving
    vψv1:: Float64 = 1.0 # right moving
    vψ1 :: Float64 = 1.0 # left moving
    # components of Az principal matrix
    az11 :: Float64
    az12 :: Float64 
    az13 :: Float64
    az21 :: Float64
    az22 :: Float64
    az23 :: Float64
    az31 :: Float64
    az32 :: Float64
    az33 :: Float64
    # source terms
    b11 :: Float64
    b12 :: Float64 
    b13 :: Float64
    b21 :: Float64
    b22 :: Float64
    b23 :: Float64
    b31 :: Float64
    b32 :: Float64
    b33 :: Float64
end

@with_kw struct noise_cibvp <: CIBVP
    noise_amp_ϕ2  :: Float64 = 0.0
    noise_amp_ψv2 :: Float64 = 0.0
    noise_amp_ψ2  :: Float64 = 0.0
    # speed of ψ2 (absolute value). To match with vψ1 (in char. frame), vψ2 = vψ1/(1+vψ1)
    vψ2 :: Float64 = 0.5 # left moving
    # components of Az principal matrix
    az11 :: Float64
    az12 :: Float64 
    az13 :: Float64
    az21 :: Float64
    az22 :: Float64
    az23 :: Float64
    az31 :: Float64
    az32 :: Float64
    az33 :: Float64
    # source terms
    b11 :: Float64
    b12 :: Float64
    b13 :: Float64
    b21 :: Float64
    b22 :: Float64
    b23 :: Float64
    b31 :: Float64
    b32 :: Float64
    b33 :: Float64
end

"""
The following functions return a real number, it is made into a grid
function in model_CCE_CCM.jl via loops
"""

# Cauchy initial data
model_CCE_CCM.ϕ1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ϕ1 * randn(T)
model_CCE_CCM.ψv1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψv1 * randn(T)
model_CCE_CCM.ψ1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψ1 * randn(T)
# Cauchy boundary data
model_CCE_CCM.ϕ1_BD(t::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ϕ1 * randn(T)
model_CCE_CCM.ψv1_BD(t::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψv1 * randn(T)
model_CCE_CCM.ψ1_BD(t::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψ1 * randn(T)
# Characteristic initial data
model_CCE_CCM.ψ2_ID(ρ::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ψ2 * randn(T)
# characteristic boundary data
model_CCE_CCM.ϕ2_BD(t::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ϕ2 * randn(T)
model_CCE_CCM.ψv2_BD(t::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ψv2 * randn(T)
model_CCE_CCM.ψ2_BD(t::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ψ2 * randn(T)

toy_model = "SYMH_SYMH_noise_t20_L2_amp"
root_dir="./run_ibvp_cibvp/"

# change D for number of points
D = 0
Nρ = (16)*2^D + 1
NX = Nρ
Nz = 16*2^D #16 coarse
# given data noise amplitude drop for fields WITHOUT derivatives in the norm:
noise_amplitude_drop_a = 0.25 
# given data noise amplitude drop for fields WITH derivatives in the norm:
noise_amplitude_drop_b = 0.125

# copy parfile in outdir
par_copy = joinpath(root_dir, toy_model, "data_$(NX)_$(Nz)/run_ibvp_cibvp_noise.jl")
mkpath(par_copy)
cp("./run_ibvp_cibvp_noise.jl", par_copy, force=true)

# parameters to be passed in the model
p = Param(
    NX = NX,
    Nρ = Nρ,
    Nz = Nz,
    ρmin = -1.0,
    ρmax = 0.0,
    Xmin = 0.0,
    Xmax = 1.0,
    tmax = 20.0, # total time of the simulation in code units
    #CFL
    cfl = 0.25,
    out_dir = joinpath(root_dir, toy_model,
                       "data_$(NX)_$(Nz)"),
    out_every = 4*2^D,
)

# the sate vector is v1 = [ϕ1, ψv1, ψ1]
ibvp = noise_ibvp(
    # for given data that converge in q-norm, a->b ϕ1 and ϕ2 amp; for H1, a->b everywhere
    noise_amp_ϕ1  = noise_amplitude_drop_a^D,
    noise_amp_ψv1 = noise_amplitude_drop_a^D,
    noise_amp_ψ1  = noise_amplitude_drop_a^D,
    # Az principal part
    az11 =  0.0,
    az12 =  1.0, # 1.0 for SYMH, 0.0 for WH
    az13 =  0.0, 
    az21 =  1.0, # 1.0 for SYMH
    az22 =  0.0,
    az23 =  0.0, 
    az31 =  0.0,
    az32 =  0.0, 
    az33 =  1.0, # 1.0 for SYMH
    # sources
    b11 =  0.0,
    b12 =  0.0,
    b13 =  0.0,
    b21 =  0.0,
    b22 =  0.0,
    b23 =  0.0,
    b31 =  0.0,
    b32 =  0.0,
    b33 =  0.0    
)

# the sate vector is v2 = [ϕ2, ψv2, ψ2]
cibvp = noise_cibvp(
    # for given data that converge in q-norm, a->b ϕ1 and ϕ2 amp; for H1, a->b everywhere
    noise_amp_ϕ2  = noise_amplitude_drop_a^D,
    noise_amp_ψv2 = noise_amplitude_drop_a^D,
    noise_amp_ψ2  = noise_amplitude_drop_a^D,
    # left-moving speed
    # vψ2 = 0.5,
    # Az principal part
    az11 =  0.0,
    az12 =  1.0, # 1.0 SYMH; 0.0 WH
    az13 =  0.0,
    az21 =  1.0,
    az22 =  0.0,
    az23 =  0.0,
    az31 =  0.0,
    az32 =  0.0,
    az33 =  1.0,
    # sources
    b11 =  0.0,
    b12 =  0.0,
    b13 =  0.0,
    b21 =  0.0,
    b22 =  0.0,
    b23 =  0.0,
    b31 =  0.0,
    b32 =  0.0,
    b33 =  0.0    
)

run_ibvp_cibvp(p, ibvp, cibvp)
