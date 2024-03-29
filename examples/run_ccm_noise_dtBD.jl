
#using model_CCE_CCM
include("../src/model_CCE_CCM_dtBD.jl")
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
#model_CCE_CCM_dtBD.
ϕ1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ϕ1 * randn(T)
#model_CCE_CCM_dtBD.
ψv1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψv1 * randn(T)
#model_CCE_CCM_dtBD.
ψ1_ID(ρ::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψ1 * randn(T)
# Cauchy boundary data
# right-moving:
#model_CCE_CCM_dtBD.
dtϕ1_BD(t::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ϕ1 * randn(T)
#model_CCE_CCM_dtBD.
dtψv1_BD(t::T, z::T, ibvp::noise_ibvp) where {T<:Real} =
    ibvp.noise_amp_ψv1 * randn(T)
# the left-moving BD for the IBVP is given by the solution to the CIBVP for CCM
# characteristic initial data
#model_CCE_CCM_dtBD.
ψ2_ID(ρ::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ψ2 * randn(T)
# characteristic boundary data
# left-moving:
#model_CCE_CCM_dtBD.
dtψ2_BD(t::T, z::T, cibvp::noise_cibvp) where {T<:Real} =
    cibvp.noise_amp_ψ2 * randn(T)
# the right-moving BD for the CIBVP is given by the solution to the IBVP

########
# PLAY #
########

# change the name according to the setup you are solving (models + given data)
toy_model = "SYMH_B1_WH_B2_noise_t20_H1_amp" 
root_dir="./run_ccm_dtBD/"

# change D for number of points
D = 4
Nρ = (16)*2^D + 1
NX = Nρ
Nz = 16*2^D #16 coarse
# given data noise amplitude drop for fields WITHOUT derivatives in the norm:
noise_amplitude_drop_a = 0.25 
# given data noise amplitude drop for fields WITH derivatives in the norm:
noise_amplitude_drop_b = 0.125

# copy parfile in outdir
par_copy = joinpath(root_dir, toy_model, "data_$(NX)_$(Nz)/run_ccm_noise_dtBD.jl")
mkpath(par_copy)
cp("./run_ccm_noise_dtBD.jl", par_copy, force=true)

# parameters to be passed in the model
p = Param(
    NX = NX,
    Nρ = Nρ,
    Nz = Nz,
    ρmin = -1.0,
    ρmax =  0.0,
    Xmin =  0.0,
    Xmax =  1.0,
    tmax =  20.0, # total time of the simulation in code units
    #CFL
    cfl = 0.25,
    out_dir = joinpath(root_dir, toy_model,
                       "data_$(NX)_$(Nz)"),
    out_every = 4*2^D,
)

# the sate vector is v1 = [ϕ1, ψv1, ψ1]
ibvp = noise_ibvp(
    # for given data that converge in q-norm, a->b ONLY below; for H1, a->b everywhere
    noise_amp_ϕ1  = noise_amplitude_drop_b^D, 
    noise_amp_ψv1 = noise_amplitude_drop_b^D,
    noise_amp_ψ1  = noise_amplitude_drop_b^D,
    # Az principal part
    az11 =  0.0,
    az12 =  1.0, # 1.0 for SYMH, 0.0 for WH
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
    b13 =  1.0, #B1
    b21 =  0.0,
    b22 =  0.0,
    b23 =  0.0,
    b31 =  0.0,
    b32 =  0.0,
    b33 =  0.0
)

# the sate vector is v2 = [ϕ2, ψv2, ψ2]
cibvp = noise_cibvp(
    # the commented out are not needed for CCM
    # noise_amp_ϕ2 = noise_amplitude_drop_a^D,
    # noise_amp_ψ2  = noise_amplitude_drop_a^D,
    noise_amp_ψ2  = noise_amplitude_drop_b^D,
    # left-moving speed
    # vψ2 = 0.5,
    # Az principal part
    az11 =  0.0,
    az12 =  0.0, # 1.0 SYMH; 0.0 WH
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
    b32 =  1.0, #B2
    b33 =  0.0  
)

run_ccm_dtBD(p, ibvp, cibvp)
