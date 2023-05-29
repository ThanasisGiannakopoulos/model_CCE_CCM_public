
#using model_CCE_CCM
include("../src/model_CCE_CCM_dtBD.jl")
using Parameters

@with_kw struct smooth_ibvp <: IBVP
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

@with_kw struct smooth_cibvp <: CIBVP
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

a1 = 1.0   #amplitude of Gaussian for IBVP
b1 = -0.5   #center of Gaussian for IBVP
c1 = 0.055 #width of Gaussian for IBVP; chosen s.t. at the ends of the domain it is less that 10^-16
function ID_gauss_IBVP(ρ::Real, z::Real)
    f =  a1*exp( -((ρ -b1)^2)/ (2.0*c1^2) ) * sin(z)
end

a2 = 1.0   #amplitude of Gaussian for IBVP
b2 = 0.5   #center of Gaussian for IBVP
c2 = 0.055 #width of Gaussian for IBVP; chosen s.t. at the ends of the domain it is less that 10^-16
function ID_gauss_CIBVP(x::Real, z::Real)
    f =  a2*exp( -((x -b2)^2)/ (2.0*c2^2) ) * sin(z)
end

function in_gauss_various_pulses(t::Real, z::Real)
    #(exp(-((t - 3)/0.5)^2) +
    # exp(-((t - 5)/0.8)^2) + exp(-((t - 8)/0.9)^2) +
    # exp(-((t - 14)/1.1)^2) + exp(-((t - 18)/1.3)^2)
    # ) * sin(z)

    #d/dt of the above =
    
    (
        - (t - 3.0)* 8.0* exp(-((t - 3.0)/0.5)^2) -
        (t - 5.0)* 3.125* exp(-((t - 5.0)/0.8)^2) -
        (t - 8.0)* 2.46914* exp(-((t - 8.0)/0.9)^2) -
        (t - 14.0)* 1.65289* exp(-((t - 14.0)/1.1)^2) -
        (t - 18.0)* 1.18343* exp(-((t - 18.0)/1.3)^2)
    ) * sin(z)

end

# Cauchy initial data
#model_CCE_CCM.
ϕ1_ID(ρ::T, z::T, ibvp::smooth_ibvp) where {T<:Real} =
    1.0*ID_gauss_IBVP(ρ, z)
#model_CCE_CCM.
ψv1_ID(ρ::T, z::T, ibvp::smooth_ibvp) where {T<:Real} =
    1.0*ID_gauss_IBVP(ρ, z)
#model_CCE_CCM.
ψ1_ID(ρ::T, z::T, ibvp::smooth_ibvp) where {T<:Real} =
    1.0*ID_gauss_IBVP(ρ, z)
# Cauchy boundary data
# right-moving:
#model_CCE_CCM.
dtϕ1_BD(t::T, z::T, ibvp::smooth_ibvp) where {T<:Real} =
    in_gauss_various_pulses(t,z)
#0.0 # it is not compatible with ID, but we neglect the diff, it is less that 10^-16
#model_CCE_CCM.
dtψv1_BD(t::T, z::T, ibvp::smooth_ibvp) where {T<:Real} =
    in_gauss_various_pulses(t,z)
    #0.0
# the left-moving BD for the IBVP is given by the solution to the CIBVP for CCM
# characteristic initial data
#model_CCE_CCM.
ψ2_ID(x::T, z::T, cibvp::smooth_cibvp) where {T<:Real} =
    1.0*ID_gauss_CIBVP(x, z)
# characteristic boundary data
# left-moving:
#model_CCE_CCM.
dtψ2_BD(t::T, z::T, cibvp::smooth_cibvp) where {T<:Real} =
    in_gauss_various_pulses(t,z)
#0.0
# the right-moving BD for the CIBVP is given by the solution to the IBVP

########
# PLAY #
########

# change the name according to the setup you are solving (models + given data)
toy_model = "SYMH_B1_WH_B2_smooth_t20" 
root_dir="./run_ccm_dtBD/"

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
par_copy = joinpath(root_dir, toy_model, "data_$(NX)_$(Nz)/run_ccm_smooth.jl")
mkpath(par_copy)
cp("./run_ccm_smooth_dtBD.jl", par_copy, force=true)

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
ibvp = smooth_ibvp(
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
cibvp = smooth_cibvp(
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
