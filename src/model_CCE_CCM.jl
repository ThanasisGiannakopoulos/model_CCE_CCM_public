module model_CCE_CCM

export Param, IBVP, CIBVP # needed to read parameters of the discretization and of the PDE structure for the IBVP and CIBVP
export run_ibvp_cibvp, run_cce, run_ccm # exports functions that are called during the run of an example

# Julia libraries needed
using Parameters
using HDF5
using Printf
using Random

@with_kw struct Param
    # discretization parameters
    NX                :: Int
    Nz                :: Int
    Nρ                :: Int
    ρmin              :: Float64 = -1.0
    ρmax              :: Float64 =  0.0
    Xmin              :: Float64 =  0.0
    Xmax              :: Float64 =  1.0
    tmax              :: Float64
    # directory to save data
    out_dir           :: String
    #CFL
    cfl               :: Float64 = 0.25
    # how often to save data
    out_every               :: Int
    compute_L2_norm_every   :: Int = 1
    compute_q_norm_every    :: Int = 1
end

"""
Cauchy fields are labeled with 1, e.g. ϕ1
Characteristic with 2, e.g. ϕ2
"""
# types used to pass parameters related to the PDE structure the IBVP and CIBVP
abstract type IBVP end
abstract type CIBVP end

"""
Functions that define the ID of ϕ1, ψv1, and ψ1 at t=0, i.e. ψ(t=0,ρ).
Needs signature ψ1_ID(ρ::T, ibvp::IBVP) where {T<:Real}, 
and similar for ϕ1 and ψv1.
"""
function ϕ1_ID end
function ψv1_ID end
function ψ1_ID end

"""
Functions that define the BCs of ϕ1, ψv1 at ρ=ρmin (right-moving)
and for ψ1 at ρ=ρmax (left-moving)
Needs signature ϕ1_BD(t::T, ibvp::IBVP) where {T<:Real}, 
and similar for ψ1 and ψv1.
"""
function ϕ1_BD end
function ψv1_BD end
function ψ1_BD end

"""
Functions that define the ID of ψ2 at u=0, i.e. ψ(u=0,x).
Needs signature ψ2_ID(r::T, cibvp::CIBVP) where {T<:Real}.
"""
function ψ2_ID end

"""
Functions that define the BCs of ϕ2, ψv2 at r=rmin (right-moving)
and for ψ2 at r=rmax (left-moving)
Needs signature ϕ2_BD(t::T, cibvp::CIBVP) where {T<:Real}, 
and similar for ψ2 and ψv2.
"""
function ϕ2_BD end
function ψv2_BD end
function ψ2_BD end

struct System{T}
    ρ    :: Vector{T}
    hρ   :: T
    X    :: Vector{T}
    hX   :: T
    z    :: Vector{T}
    hz   :: T
end
function System(p::Param)
    # the -1 point is to have a compact domain; that includes the max value
    hρ  = (p.ρmax - p.ρmin) / (p.Nρ-1)
    ρ   = range(p.ρmin, p.ρmax, length=p.Nρ)
    hX  = (p.Xmax - p.Xmin) / (p.NX-1)
    X   = range(p.Xmin, p.Xmax, length=p.NX)
    hz  = 2.0*π / (p.Nz) # remove last point; 0=2π (periodic)
    z   = range(0, 2.0*π-hz, length=p.Nz)
    System{typeof(hX)}(ρ, hρ, X, hX, z, hz)
end

# load the discretized derivative operators
include("./operators.jl")

# functions for initial data
# IBVP
function initialize_ϕ1(sys::System{T}, ibvp::IBVP) where {T}
    Nρ = length(sys.ρ)
    Nz = length(sys.z)

    f0 = zeros(T, Nρ, Nz)

    @inbounds for j in 1:Nz
        @inbounds for i in 1:Nρ
            f0[i,j] = ϕ1_ID(sys.ρ[i], sys.z[j], ibvp)
        end
    end
    f0
end
function initialize_ψv1(sys::System{T}, ibvp::IBVP) where {T}
    Nρ = length(sys.ρ)
    Nz = length(sys.z)

    f0 = zeros(T, Nρ, Nz)

    @inbounds for j in 1:Nz
        @inbounds for i in 1:Nρ
            f0[i,j] = ψv1_ID(sys.ρ[i], sys.z[j], ibvp)
        end
    end
    f0
end
function initialize_ψ1(sys::System{T}, ibvp::IBVP) where {T}
    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    
    f0 = zeros(T, Nρ, Nz)

    @inbounds for j in 1:Nz
        @inbounds for i in 1:Nρ
            f0[i,j] = ψ1_ID(sys.ρ[i], sys.z[j], ibvp)
        end
    end
    f0
end

# CIBVP
function initialize_ψ2(sys::System{T}, cibvp::CIBVP) where {T}
    NX = length(sys.X)
    Nz = length(sys.z)

    f0 = zeros(T, NX, Nz)

    @inbounds for j in 1:Nz
        @inbounds for i in 1:NX
            f0[i,j] = ψ2_ID(sys.X[i], sys.z[j], cibvp)
        end
    end
    f0
end

# function to write data
function write_2D(it::Int, t::Float64, data_dir::String,
                  ϕ1::Matrix, ψv1::Matrix, ψ1::Matrix,
                  ϕ2::Matrix, ψv2::Matrix, ψ2::Matrix)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "data_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "ϕ1" ,  ϕ1)
        write(file, "ψv1",  ψv1)
        write(file, "ψ1" ,  ψ1)
        write(file, "ϕ2" ,  ϕ2)
        write(file, "ψv2",  ψv2)
        write(file, "ψ2" ,  ψ2)
        attrs(file)["time"] = t
    end
    nothing
end

"""
intrinsic equations of the CIBVP
v[:,1] = ϕ_2(x_c, z), 
and
v[:,2] = ψ_v2(x_c, z), 
for some fixed x_c.
dv is the derivative ∂_x*(ϕ_2, ψ_v2)^T.

cibvp.aij and cibvp.bij, are parameters that control the structure of
the principal part and sources (respectively) of the intrisic
equations for the CIBVP.
 """

# mutates dv
function intrinsic_eq_rhs!(dv::Matrix, v::Matrix,
                           S::Vector,
                           sys::System, cibvp::CIBVP)

    dv[:,1] .= cibvp.az11*Dz(v[:,1],sys) .+ cibvp.az12*Dz(v[:,2],sys) .+
        cibvp.az13*Dz(S[:],sys) .+
        cibvp.b11*v[:,1] .+ cibvp.b12*v[:,2] .+ cibvp.b13*S[:]
    dv[:,2] .= cibvp.az21*Dz(v[:,1],sys) .+ cibvp.az22*Dz(v[:,2],sys) .+
        cibvp.az23*Dz(S[:],sys) .+
        cibvp.b21*v[:,1] .+ cibvp.b22*v[:,2] .+ cibvp.b23*S[:]

end

"""
trapezoidal rule to solve the explicitly solve the instrinsic
equations of the CIBVP
"""

# mutates vf
function trapez!(hX::Float64,
                 vf::Matrix,
                 v0::Matrix, v1::Matrix, v2::Matrix,
                 S0::Vector, # source at x_i-1
                 S1::Vector, # source at x_i
                 sys::System, cibvp::CIBVP)

    intrinsic_eq_rhs!(v1, v0, S0, sys, cibvp)
    intrinsic_eq_rhs!(v2, v0 .+ hX*v1, S1, sys, cibvp)
    
    vf .= v0 .+ 0.5*hX*(v1 .+ v2)
    #nothing
    vf
end

"""
The boundary data treatment is different when the IBVP and CIBVP are solved seperately, in CCE and in CCM.
"""

##################################################################################################
# IBVP and CIBVP seperately
##################################################################################################

# mutates dv and null_v
function setup_rhs_ibvp_cibvp!(t::Float64,
                    dv::Array{Float64, 3}, # dv/dt for  [ϕ1, ψv1, ψ1, ψ2]; d/du for characteristic
                    v::Array{Float64, 3}, # vector v = [ϕ1, ψv1, ψ1, ψ2]
                    null_v::Array{Float64, 3}, # vector null_v = [ϕ2, ψv2]
                    null_v1::Matrix, null_v2::Matrix, # used for null integration
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # boundary data (BD)
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t, sys.z[j], ibvp)
        # incoming boundary data from the right boundary of the characteristic
        v[end,j,4] = ψ2_BD(t, sys.z[j], cibvp)

        """ This part is different for CCM, CCE, and IBVP_CIBVP"""
        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(t, sys.z[j], ibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        
        # BD for null integration
        null_v[1,j,1] =  ϕ2_BD(t, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        null_v[1,j,2] =  ψv2_BD(t, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP
    end
    
    # null integration
    # updates null_v[:,:,1] and null_v[:,:,2]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    #dϕ1/dt; right moving
    dv[:,:,1] .= -ibvp.vϕ1*Dρ(v[:,:,1], sys) .+
        ibvp.az11*Dz(v[:,:,1], sys) .+ibvp.az12*Dz(v[:,:,2], sys) .+ibvp.az13*Dz(v[:,:,3], sys) .+
        ibvp.b11*v[:,:,1] .+ ibvp.b12*v[:,:,2] .+ ibvp.b13*v[:,:,3] 
    #dψv1/dt; right moving
    dv[:,:,2] .= -ibvp.vψv1*Dρ(v[:,:,2], sys) .+
        ibvp.az21*Dz(v[:,:,1], sys) .+ibvp.az22*Dz(v[:,:,2], sys) .+ibvp.az23*Dz(v[:,:,3], sys) .+
        ibvp.b21*v[:,:,1] .+ ibvp.b22*v[:,:,2] .+ ibvp.b23*v[:,:,3] 
    #dψ1/dt; left moving
    dv[:,:,3] .= ibvp.vψ1*Dρ(v[:,:,3], sys) .+
        ibvp.az31*Dz(v[:,:,1], sys) .+ibvp.az32*Dz(v[:,:,2], sys) .+ibvp.az33*Dz(v[:,:,3], sys) .+
        ibvp.b31*v[:,:,1] .+ ibvp.b32*v[:,:,2] .+ ibvp.b33*v[:,:,3]
    #dψ2/du in characteristic; left moving
    dv[:,:,4] .= cibvp.vψ2* DX(v[:,:,4], sys) .+
        cibvp.az13*Dz(null_v[:,:,1], sys) .+ cibvp.az23*Dz(null_v[:,:,2], sys) .+
        cibvp.az33*Dz(v[:,:,4], sys) .+
        cibvp.b31*null_v[:,:,1] .+ cibvp.b32*null_v[:,:,2] .+ cibvp.b33*v[:,:,4]

    nothing
end

# the state vector and its derivative is a 3-tensor for a function f(x,z,t)
# with type Array{Float64, 3}
# mutates v
# the time evolution is performed with Runge-Kutta 4 (below)
function time_evol_ibvp_cibvp!(t::Float64, dt::Float64,
                    v::Array{Float64, 3},
                    v1::Array{Float64, 3}, v2::Array{Float64, 3},
                    v3::Array{Float64, 3}, v4::Array{Float64, 3},
                    null_v::Array{Float64, 3},
                    null_v1::Matrix, null_v2::Matrix,
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX
    # time evol with RK4:
    #v1 .=
    setup_rhs_ibvp_cibvp!(t, v1, v,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v2 .=
    setup_rhs_ibvp_cibvp!(t, v2, v .+ 0.5.*dt.*v1,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v3 .=
    setup_rhs_ibvp_cibvp!(t, v3, v .+ 0.5.*dt.*v2,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v4 .=
    setup_rhs_ibvp_cibvp!(t, v4, v .+ dt.*v3,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)

    # time evolve; the standard RK4 scheme
    v .= v .+ dt*(v1 .+ 2.0.*v2 .+ 2.0.*v3 .+ v4)/6.0

    # overwrite TRUE boundary data; these data will be saved
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t+dt, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t+dt, sys.z[j], ibvp)
        ## BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(t+dt, sys.z[j], ibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        # incoming boundary data from future null infinity
        v[end,j,4] = ψ2_BD(t+dt, sys.z[j], cibvp)
    
        # update ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = ϕ2_BD(t+dt, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        null_v[1,j,2] = ψv2_BD(t+dt, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP

    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    nothing
end

# function that performs the time evolution
function run_ibvp_cibvp(p::Param, ibvp::IBVP, cibvp::CIBVP)

    # pass the parameters of the system
    sys = System(p)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # create the folders where data are saved
    data_dir = mkpath(p.out_dir)

    # timestep
    dt = p.cfl* minimum([sys.hX, sys.hρ])
    tmax = p.tmax

    # write the ID as a vector for the coupled PDE system
    # state vector v = [ϕ1, ψv1, ψ1, ψ2]
    v  = zeros(p.Nρ, p.Nz, 4)
    v1 = zeros(p.Nρ, p.Nz, 4)
    v2 = zeros(p.Nρ, p.Nz, 4)
    v3 = zeros(p.Nρ, p.Nz, 4)
    v4 = zeros(p.Nρ, p.Nz, 4)
    # state vector null_v = [ϕ2, ψv2]
    null_v = zeros(p.NX, p.Nz, 2)
    # for null integration
    null_v1 = zeros(p.Nz,2)
    null_v2 = zeros(p.Nz,2)

    # set initial data 
    v[:,:,1] .= initialize_ϕ1(sys, ibvp) # ϕ1
    v[:,:,2] .= initialize_ψv1(sys, ibvp) # ψv1
    v[:,:,3] .= initialize_ψ1(sys, ibvp) # ψ1
    v[:,:,4] .= initialize_ψ2(sys, cibvp) # ψ2

    # boundary data
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(0.0, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(0.0, sys.z[j], ibvp)
        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(0.0, sys.z[j], ibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        # incoming boundary data from future null infinity
        v[end,j,4] = ψ2_BD(0.0, sys.z[j], cibvp)
    
        # find initial ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = ϕ2_BD(0.0, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        null_v[1,j,2] = ψv2_BD(0.0, sys.z[j], cibvp) # prescribed BD from example file; it does NOT come from the CIBVP
        
    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:p.NX
        null_v[i,:,:] .= trapez!(sys.hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    # write the coordinates. if there's a problem with this for windows users,
    # it can be done differently (with joinpath, or similar)
    h5write(data_dir*"/ρ.h5", "ρ", sys.ρ)
    h5write(data_dir*"/X.h5", "X", sys.X)
    h5write(data_dir*"/z.h5", "z", sys.z)

    it = 0
    t  = 0.0
    
    # save initial data
    write_2D(it, t, data_dir,
             v[:,:,1], v[:,:,2], v[:,:,3],
             null_v[:,:,1], null_v[:,:,2], v[:,:,4])
    # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
    
    #println("number of threads = ", Threads.nthreads())
    println("Running IBVP and CIBVP separately")
    println("-------------------------------------------------------------")
    println("Iteration      Time |           ψ2            |           ψ1           |           ϕ2            ")
    println("                    |    minimum      maximum |    minimum      maximum |    minimum      maximum |")
    println("-------------------------------------------------------------")
    
    @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])
    
    # start time evolution
    while t<tmax
        it += 1
        t  += dt
        
        time_evol_ibvp_cibvp!(t, dt,
                   v, v1, v2, v3, v4,
                   null_v,
                   null_v1, null_v2,
                   sys, ibvp, cibvp)
        
        @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])

        # save data
        if it % p.out_every == 0
            write_2D(it, t, data_dir,
                     v[:,:,1], v[:,:,2], v[:,:,3],
                     null_v[:,:,1], null_v[:,:,2], v[:,:,4])
            # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
        end
        
    end # end evolution

    println("-------------------------------------------------------------")
    println("Done.")
end # end run_ibvp_cibvp

##################################################################################################
# CCE
##################################################################################################

function setup_rhs_cce!(t::Float64,
                    dv::Array{Float64, 3}, # dv/dt for  [ϕ1, ψv1, ψ1, ψ2]; d/du for characteristic
                    v::Array{Float64, 3}, # vector v = [ϕ1, ψv1, ψ1, ψ2]
                    null_v::Array{Float64, 3}, # vector null_v = [ϕ2, ψv2]
                    null_v1::Matrix, null_v2::Matrix, # used for null integration
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # boundary data (BD)
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t, sys.z[j], ibvp)
        # incoming boundary data from the right boundary of the characteristic domain
        v[end,j,4] = ψ2_BD(t, sys.z[j], cibvp)

        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(t, sys.z[j], ibvp) # for CCE: prescribed BD for the left-moving field of the IBVP
        
        # BD for null integration
        null_v[1,j,1] =  v[end,j,1] # for CCE: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] =  v[end,j,2] # for CCE: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
    end
    
    # null integration
    # updates null_v[:,:,1] and null_v[:,:,2]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    #dϕ1/dt; right moving
    dv[:,:,1] .= -ibvp.vϕ1*Dρ(v[:,:,1], sys) .+
        ibvp.az11*Dz(v[:,:,1], sys) .+ibvp.az12*Dz(v[:,:,2], sys) .+ibvp.az13*Dz(v[:,:,3], sys) .+
        ibvp.b11*v[:,:,1] .+ ibvp.b12*v[:,:,2] .+ ibvp.b13*v[:,:,3] 
    #dψv1/dt; right moving
    dv[:,:,2] .= -ibvp.vψv1*Dρ(v[:,:,2], sys) .+
        ibvp.az21*Dz(v[:,:,1], sys) .+ibvp.az22*Dz(v[:,:,2], sys) .+ibvp.az23*Dz(v[:,:,3], sys) .+
        ibvp.b21*v[:,:,1] .+ ibvp.b22*v[:,:,2] .+ ibvp.b23*v[:,:,3] 
    #dψ1/dt; left moving
    dv[:,:,3] .= ibvp.vψ1*Dρ(v[:,:,3], sys) .+
        ibvp.az31*Dz(v[:,:,1], sys) .+ibvp.az32*Dz(v[:,:,2], sys) .+ibvp.az33*Dz(v[:,:,3], sys) .+
        ibvp.b31*v[:,:,1] .+ ibvp.b32*v[:,:,2] .+ ibvp.b33*v[:,:,3]
    #dψ2/du in characteristic; left moving
    dv[:,:,4] .= cibvp.vψ2* DX(v[:,:,4], sys) .+
        cibvp.az13*Dz(null_v[:,:,1], sys) .+ cibvp.az23*Dz(null_v[:,:,2], sys) .+
        cibvp.az33*Dz(v[:,:,4], sys) .+
        cibvp.b31*null_v[:,:,1] .+ cibvp.b32*null_v[:,:,2] .+ cibvp.b33*v[:,:,4]

    nothing
end

function time_evol_cce!(t::Float64, dt::Float64,
                    v::Array{Float64, 3},
                    v1::Array{Float64, 3}, v2::Array{Float64, 3},
                    v3::Array{Float64, 3}, v4::Array{Float64, 3},
                    null_v::Array{Float64, 3},
                    null_v1::Matrix, null_v2::Matrix,
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX
    # time evol with RK4:
    #v1 .=
    setup_rhs_cce!(t, v1, v,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v2 .=
    setup_rhs_cce!(t, v2, v .+ 0.5.*dt.*v1,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v3 .=
    setup_rhs_cce!(t, v3, v .+ 0.5.*dt.*v2,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v4 .=
    setup_rhs_cce!(t, v4, v .+ dt.*v3,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)

    # time evolve; the standard RK4 scheme
    v .= v .+ dt*(v1 .+ 2.0.*v2 .+ 2.0.*v3 .+ v4)/6.0

    # overwrite TRUE boundary data; these data will be saved
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t+dt, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t+dt, sys.z[j], ibvp)
        ## BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(t+dt, sys.z[j], ibvp) # for CCE: prescribed BD for the left-moving field of the IBVP
        # incoming boundary data from the right boundary of the CIBVP
        v[end,j,4] = ψ2_BD(t+dt, sys.z[j], cibvp)
    
        # update ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = v[end,j,1] # for CCE: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] = v[end,j,2] # for CCE: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain

    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    nothing
end

function run_cce(p::Param, ibvp::IBVP, cibvp::CIBVP)

    # pass the parameters of the system
    sys = System(p)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # create the folders where data are saved
    data_dir = mkpath(p.out_dir)

    # timestep
    dt = p.cfl* minimum([sys.hX, sys.hρ])
    tmax = p.tmax

    # write the ID as a vector for the coupled PDE system
    # state vector v = [ϕ1, ψv1, ψ1, ψ2]
    v  = zeros(p.Nρ, p.Nz, 4)
    v1 = zeros(p.Nρ, p.Nz, 4)
    v2 = zeros(p.Nρ, p.Nz, 4)
    v3 = zeros(p.Nρ, p.Nz, 4)
    v4 = zeros(p.Nρ, p.Nz, 4)
    # state vector null_v = [ϕ2, ψv2]
    null_v = zeros(p.NX, p.Nz, 2)
    # for null integration
    null_v1 = zeros(p.Nz,2)
    null_v2 = zeros(p.Nz,2)

    # set initial data 
    v[:,:,1] .= initialize_ϕ1(sys, ibvp) # ϕ1
    v[:,:,2] .= initialize_ψv1(sys, ibvp) # ψv1
    v[:,:,3] .= initialize_ψ1(sys, ibvp) # ψ1
    v[:,:,4] .= initialize_ψ2(sys, cibvp) # ψ2

    # boundary data
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(0.0, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(0.0, sys.z[j], ibvp)
        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = ψ1_BD(0.0, sys.z[j], ibvp) # for CCE: prescribed BD for the left-moving field of the IBVP
        # incoming boundary data from the right boundary of the CIBVP
        v[end,j,4] = ψ2_BD(0.0, sys.z[j], cibvp)
    
        # find initial ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = v[end,j,1] # for CCE: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] = v[end,j,2] # for CCE: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        
    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:p.NX
        null_v[i,:,:] .= trapez!(sys.hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    # write the coordinates. if there's a problem with this for windows users,
    # it can be done differently (with joinpath, or similar)
    h5write(data_dir*"/ρ.h5", "ρ", sys.ρ)
    h5write(data_dir*"/X.h5", "X", sys.X)
    h5write(data_dir*"/z.h5", "z", sys.z)

    it = 0
    t  = 0.0
    
    # save initial data
    write_2D(it, t, data_dir,
             v[:,:,1], v[:,:,2], v[:,:,3],
             null_v[:,:,1], null_v[:,:,2], v[:,:,4])
    # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
    
    #println("number of threads = ", Threads.nthreads())
    println("Running CCE")
        println("-------------------------------------------------------------")
    println("Iteration      Time |           ψ2            |           ψ1           |           ϕ2            ")
    println("                    |    minimum      maximum |    minimum      maximum |    minimum      maximum |")
    println("-------------------------------------------------------------")
    
    @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])
    
    # start time evolution
    while t<tmax
        it += 1
        t  += dt
        
        time_evol_cce!(t, dt,
                   v, v1, v2, v3, v4,
                   null_v,
                   null_v1, null_v2,
                   sys, ibvp, cibvp)
        
        @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])

        # save data
        if it % p.out_every == 0
            write_2D(it, t, data_dir,
                     v[:,:,1], v[:,:,2], v[:,:,3],
                     null_v[:,:,1], null_v[:,:,2], v[:,:,4])
            # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
        end
        
    end # end evolution

    println("-------------------------------------------------------------")
    println("Done.")
end # end run_cce

##################################################################################################
# CCM
##################################################################################################

function setup_rhs_ccm!(t::Float64,
                    dv::Array{Float64, 3}, # dv/dt for  [ϕ1, ψv1, ψ1, ψ2]; d/du for characteristic
                    v::Array{Float64, 3}, # vector v = [ϕ1, ψv1, ψ1, ψ2]
                    null_v::Array{Float64, 3}, # vector null_v = [ϕ2, ψv2]
                    null_v1::Matrix, null_v2::Matrix, # used for null integration
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # boundary data (BD)
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t, sys.z[j], ibvp)
        # incoming boundary data from the right boundary of the characteristic
        v[end,j,4] = ψ2_BD(t, sys.z[j], cibvp)

        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = v[1,j,4] # for CCM: ψ1[end,j] = ψ2[1,j]; left-moving, pass information from the characteristic to the Cauchy domain
        
        # BD for null integration
        null_v[1,j,1] =  v[end,j,1] # for CCM: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] =  v[end,j,2] # for CCM: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
    end
    
    # null integration
    # updates null_v[:,:,1] and null_v[:,:,2]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    #dϕ1/dt; right moving
    dv[:,:,1] .= -ibvp.vϕ1*Dρ(v[:,:,1], sys) .+
        ibvp.az11*Dz(v[:,:,1], sys) .+ibvp.az12*Dz(v[:,:,2], sys) .+ibvp.az13*Dz(v[:,:,3], sys) .+
        ibvp.b11*v[:,:,1] .+ ibvp.b12*v[:,:,2] .+ ibvp.b13*v[:,:,3] 
    #dψv1/dt; right moving
    dv[:,:,2] .= -ibvp.vψv1*Dρ(v[:,:,2], sys) .+
        ibvp.az21*Dz(v[:,:,1], sys) .+ibvp.az22*Dz(v[:,:,2], sys) .+ibvp.az23*Dz(v[:,:,3], sys) .+
        ibvp.b21*v[:,:,1] .+ ibvp.b22*v[:,:,2] .+ ibvp.b23*v[:,:,3] 
    #dψ1/dt; left moving
    dv[:,:,3] .= ibvp.vψ1*Dρ(v[:,:,3], sys) .+
        ibvp.az31*Dz(v[:,:,1], sys) .+ibvp.az32*Dz(v[:,:,2], sys) .+ibvp.az33*Dz(v[:,:,3], sys) .+
        ibvp.b31*v[:,:,1] .+ ibvp.b32*v[:,:,2] .+ ibvp.b33*v[:,:,3]
    #dψ2/du in characteristic; left moving
    dv[:,:,4] .= cibvp.vψ2* DX(v[:,:,4], sys) .+
        cibvp.az13*Dz(null_v[:,:,1], sys) .+ cibvp.az23*Dz(null_v[:,:,2], sys) .+
        cibvp.az33*Dz(v[:,:,4], sys) .+
        cibvp.b31*null_v[:,:,1] .+ cibvp.b32*null_v[:,:,2] .+ cibvp.b33*v[:,:,4]

    nothing
end

function time_evol_ccm!(t::Float64, dt::Float64,
                    v::Array{Float64, 3},
                    v1::Array{Float64, 3}, v2::Array{Float64, 3},
                    v3::Array{Float64, 3}, v4::Array{Float64, 3},
                    null_v::Array{Float64, 3},
                    null_v1::Matrix, null_v2::Matrix,
                    sys::System, ibvp::IBVP, cibvp::CIBVP)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX
    # time evol with RK4:
    #v1 .=
    setup_rhs_ccm!(t, v1, v,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v2 .=
    setup_rhs_ccm!(t, v2, v .+ 0.5.*dt.*v1,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v3 .=
    setup_rhs_ccm!(t, v3, v .+ 0.5.*dt.*v2,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)
    #v4 .=
    setup_rhs_ccm!(t, v4, v .+ dt.*v3,
               null_v,
               null_v1, null_v2,
               sys, ibvp, cibvp)

    # time evolve; the standard RK4 scheme
    v .= v .+ dt*(v1 .+ 2.0.*v2 .+ 2.0.*v3 .+ v4)/6.0

    # overwrite TRUE boundary data; these data will be saved
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(t+dt, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(t+dt, sys.z[j], ibvp)
        ## BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = v[1,j,4] # for CCM: ψ1[end,j] = ψ2[1,j]; left-moving, pass information from the characteristic to the Cauchy domain
        # incoming boundary data from future null infinity
        v[end,j,4] = ψ2_BD(t+dt, sys.z[j], cibvp)
    
        # update ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = v[end,j,1] # for CCM: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] = v[end,j,2] # for CCM: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain

    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:NX
        null_v[i,:,:] .= trapez!(hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    nothing
end

# function that performs the time evolution
function run_ccm(p::Param, ibvp::IBVP, cibvp::CIBVP)

    # pass the parameters of the system
    sys = System(p)

    Nρ = length(sys.ρ)
    Nz = length(sys.z)
    NX = Nρ
    hX = sys.hX

    # create the folders where data are saved
    data_dir = mkpath(p.out_dir)

    # timestep
    dt = p.cfl* minimum([sys.hX, sys.hρ])
    tmax = p.tmax

    # write the ID as a vector for the coupled PDE system
    # state vector v = [ϕ1, ψv1, ψ1, ψ2]
    v  = zeros(p.Nρ, p.Nz, 4)
    v1 = zeros(p.Nρ, p.Nz, 4)
    v2 = zeros(p.Nρ, p.Nz, 4)
    v3 = zeros(p.Nρ, p.Nz, 4)
    v4 = zeros(p.Nρ, p.Nz, 4)
    # state vector null_v = [ϕ2, ψv2]
    null_v = zeros(p.NX, p.Nz, 2)
    # for null integration
    null_v1 = zeros(p.Nz,2)
    null_v2 = zeros(p.Nz,2)

    # set initial data 
    v[:,:,1] .= initialize_ϕ1(sys, ibvp) # ϕ1
    v[:,:,2] .= initialize_ψv1(sys, ibvp) # ψv1
    v[:,:,3] .= initialize_ψ1(sys, ibvp) # ψ1
    v[:,:,4] .= initialize_ψ2(sys, cibvp) # ψ2

    # boundary data
    @inbounds for j in 1:Nz
        # BD from left boundary of Cauchy
        v[1,j,1] = ϕ1_BD(0.0, sys.z[j], ibvp)
        v[1,j,2] = ψv1_BD(0.0, sys.z[j], ibvp)
        # BD from right boundary of Cauchy, that touches characteristic domain
        v[end,j,3] = v[1,j,4] # for CCM: ψ1[end,j] = ψ2[1,j]; left-moving, pass information from the characteristic to the Cauchy domain
        # incoming boundary data from future null infinity
        v[end,j,4] = ψ2_BD(0.0, sys.z[j], cibvp)
    
        # find initial ϕ2 and ψv2 as well
        # BD for null integration
        null_v[1,j,1] = v[end,j,1] # for CCM: ϕ2[1,:] = ϕ1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
        null_v[1,j,2] = v[end,j,2] # for CCM: ψv2[1,:] = ψv1[end,:]; right-moving, pass information from the Cauchy to the characteristic domain
    end
    
    # null integration
    # updates v[:,:,5] and v[:,:,6]
    @inbounds for i in 2:p.NX
        null_v[i,:,:] .= trapez!(sys.hX,
                                 null_v[i,:,:], # vf
                                 null_v[i-1,:,:], # v0
                                 null_v1, null_v2,
                                 v[i-1,:,4], # ψ2 source at x_i-1
                                 v[i,:,4], # ψ2 source at x_i
                                 sys, cibvp)
    end
    
    # write the coordinates. if there's a problem with this for windows users,
    # it can be done differently (with joinpath, or similar)
    h5write(data_dir*"/ρ.h5", "ρ", sys.ρ)
    h5write(data_dir*"/X.h5", "X", sys.X)
    h5write(data_dir*"/z.h5", "z", sys.z)

    it = 0
    t  = 0.0
    
    # save initial data
    write_2D(it, t, data_dir,
             v[:,:,1], v[:,:,2], v[:,:,3],
             null_v[:,:,1], null_v[:,:,2], v[:,:,4])
    # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
    
    #println("number of threads = ", Threads.nthreads())
    println("Running CCM")
    println("-------------------------------------------------------------")
    println("Iteration      Time |           ψ2            |           ψ1           |           ϕ2            ")
    println("                    |    minimum      maximum |    minimum      maximum |    minimum      maximum |")
    println("-------------------------------------------------------------")
    
    @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])
    
    # start time evolution
    while t<tmax
        it += 1
        t  += dt
        
        time_evol_ccm!(t, dt,
                   v, v1, v2, v3, v4,
                   null_v,
                   null_v1, null_v2,
                   sys, ibvp, cibvp)
        
        @printf "%9d %9.3f |  %9.4g    %9.4g|  %9.4g    %9.4g|  %9.4g    %9.4g\n" it t minimum(v[:,:,4]) maximum(v[:,:,4]) minimum(v[:,:,3]) maximum(v[:,:,3]) minimum(null_v[:,:,1]) maximum(null_v[:,:,1])

        # save data
        if it % p.out_every == 0
            write_2D(it, t, data_dir,
                     v[:,:,1], v[:,:,2], v[:,:,3],
                     null_v[:,:,1], null_v[:,:,2], v[:,:,4])
            # v = [ϕ1, ψv1, ψ1, ψ2, ψv2, ϕ2]
        end
        
    end # end evolution

    println("-------------------------------------------------------------")
    println("Done.")
end # end run_ccm

end # module
