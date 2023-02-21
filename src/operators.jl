
""" 
2nd order accurate finite differences for 1st order derivatives
along the different directions, in matrix form 
"""

# z direction (periodic)
function Dz!(f_z::Matrix, f::Matrix, sys::System)
    NX, Nz = size(f)
    odz2 = 0.5 / sys.hz

    @inbounds for j in 2:Nz-1
        @inbounds for i in 1:NX
            f_z[i,j] = (f[i,j+1] - f[i,j-1]) * odz2
        end
    end

    @inbounds for i in 1:NX
        f_z[i,1]   = (f[i,2] - f[i,end]) * odz2
        f_z[i,end] = (f[i,1] - f[i,end-1]) * odz2
    end

    f_z
end
# now acting on vectors
function Dz!(f_z::Vector, f::Vector, sys::System)
    Nz = length(f)
    odz2 = 0.5 / sys.hz

    @inbounds for j in 2:Nz-1
        f_z[j] = (f[j+1] - f[j-1]) * odz2
    end

    f_z[1]   = (f[2] - f[end]) * odz2
    f_z[end] = (f[1] - f[end-1]) * odz2

    f_z
end
function Dz(f, sys::System)
    f_z = similar(f)
    Dz!(f_z, f, sys)
end

"""
2nd order accurate finite difference operator for 1st order
derivatives along a non-periodic direction. Forward and backward 2nd
order accurate finite differences are used on the first and last
points of the grid respectively
"""

function DX!(f_X::Matrix,
 f::Matrix, sys::System) NX, Nz = size(f) odX2 = 0.5 / sys.hX

    @inbounds for j in 1:Nz
        @inbounds for i in 2:NX-1
            f_X[i,j] = (f[i+1,j] - f[i-1,j]) * odX2
        end
    end

    @inbounds for j in 1:Nz
        #truncation error matched
        f_X[1,j]   = ( -4.0* f[1,j] + 7.0*f[2,j] - 4.0* f[3,j] + f[4,j]) * odX2
        f_X[end,j] = ( 4.0* f[end,j] - 7.0*f[end-1,j] + 4.0* f[end-2,j] - f[end-3,j]) * odX2
    end

    f_X
end
function DX(f, sys::System)
    f_X = similar(f)
    DX!(f_X, f, sys)
end

# ρ direction
function  Dρ!(f_ρ::Matrix, f::Matrix, sys::System)
    Nρ, Nz = size(f)
    odρ2 = 0.5 / sys.hρ

    @inbounds for j in 1:Nz
        @inbounds for i in 2:Nρ-1
            f_ρ[i,j] = (f[i+1,j] - f[i-1,j]) * odρ2
        end
    end
    
    @inbounds for j in 1:Nz
        #trunctation error matched
        f_ρ[1,j]  = (-4.0* f[1,j] + 7.0*f[2,j] - 4.0* f[3,j] + f[4,j]) * odρ2
        f_ρ[end,j] = (4.0* f[end,j] - 7.0*f[end-1,j] + 4.0* f[end-2,j] - f[end-3,j]) * odρ2 
    end
    
    f_ρ
end
function Dρ(f, sys::System)
    f_ρ = similar(f)
    Dρ!(f_ρ, f, sys)
end
