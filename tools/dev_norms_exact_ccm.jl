
using HDF5
using DelimitedFiles
using Printf
using ProgressBars

# derivatives
# radial; same for X and ρ
function DX!(f_X::Matrix, f::Matrix, hX::Float64)
    NX, Nz = size(f)
    odX2 = 0.5 / hX

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
function DX(f, hX::Float64)
    f_X = similar(f)
    DX!(f_X, f, hX)
end
# angular
function Dz!(f_z::Matrix, f::Matrix, hz::Float64)
    NX, Nz = size(f)
    odz2 = 0.5 / hz

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
function Dz(f, hz::Float64)
    f_z = similar(f)
    Dz!(f_z, f, hz)
end

function list_h5_files(foldername::String; prefix::String="data_")
    path     = abspath(foldername)
    allfiles = readdir(path)

    Ns = length(prefix)

    its_names = Tuple[]
    # append only the files whose names start with the given prefix
    for file in allfiles
        try
            if (file[1:Ns] == prefix && file[end-2:end] == ".h5")
                # extract iteration
                it_str = file[Ns+1:end-3]
                fullname = joinpath(path, file)
                # add to list of tuples with iteration and name
                push!(its_names, (parse(Int64, it_str), fullname))
            end
        catch ex
            if isa(ex, BoundsError)
                # probably triggered by string comparison; do nothing
            else
                throw(ex)
            end
        end
    end

    # sort according to iteration
    sort!(its_names)
    # and extract the list of filenames and iterations
    filenames = [name for (it, name) in its_names]
    its       = [it for (it, name) in its_names]

    (its, filenames)
end

#dt0 is needed for calculation on common timestep with coarsest resol only
function dev_t_func(dir, dt0)
    x = h5read(dir * "/X.h5", "X")
    ρ = h5read(dir * "/ρ.h5", "ρ")
    z = h5read(dir * "/z.h5", "z")
    Nρ = length(ρ)
    Nx = length(x)
    Nz = length(z)
    
    # grid spacing
    dx = x[2] - x[1] 
    dρ = ρ[2] - ρ[1] 
    dz = z[2] - z[1] 
    
    # list all available iterations (and corresponding files)
    (its, all_filenames) = list_h5_files(dir, prefix="data_")
    
    filenames = all_filenames[:]
    
    Nf           = length(filenames)
    tt           = zeros(Nf)

    # full
    q  = zeros(Nf)
    H1 = zeros(Nf)
    
    # initiate the grid function for the outgoing norm
    q_out_x = zeros(length(x))
    h1_out_x = zeros(length(x))

    for it in ProgressBar(1:Nf)
        file = filenames[it]

        # cauchy
        ϕ1  = h5read(file, "ϕ1")
        ψv1 = h5read(file, "ψv1")
        ψ1  = h5read(file, "ψ1")
        # radial derivatives
        Dρϕ1  = DX(ϕ1, dρ) 
        Dρψv1 = DX(ψv1, dρ)
        Dρψ1  = DX(ψ1, dρ)
        # angular derivatives
        Dzϕ1  = Dz(ϕ1, dz) 
        Dzψv1 = Dz(ψv1, dz)
        Dzψ1  = Dz(ψ1, dz)

        # charteristic
        ϕ2  = h5read(file, "ϕ2")
        ψv2 = h5read(file, "ψv2")
        ψ2  = h5read(file, "ψ2")
        # null derivatives
        Dxϕ2  = DX(ϕ2, dx) 
        Dxψv2 = DX(ψv2, dx)
        Dxψ2  = DX(ψ2, dx)
        # angular derivatives
        Dzϕ2  = Dz(ϕ2, dz) 
        Dzψv2 = Dz(ψv2, dz)
        Dzψ2  = Dz(ψ2, dz)

        # time
        t  = h5readattr(file, "/")["time"]
        
        # time saved
        tt[it]  = t

        # get the characteristi outgoing part of the norm
        q_out_x += dt0*dz*sum( ψv2.*ψv2 .+ ϕ2.*ϕ2 .+ Dzϕ2.*Dzϕ2, dims=2)
        h1_out_x += dt0*dz*sum( ψv2.*ψv2 .+ ϕ2.*ϕ2 .+
                                Dxψv2.*Dxψv2 .+ Dxϕ2.*Dxϕ2 .+
                                Dzψv2.*Dzψv2 .+ Dzϕ2.*Dzϕ2, dims=2)
        

        # q norm
        q[it] = dx*dz*sum(ψ2.*ψ2) + maximum(q_out_x) +
            dρ*dz*sum(ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1) +
            dt0*dz*sum(ψ1[1,:].*ψ1[1,:])

        # H1 norm
        H1[it] = dx*dz*sum(ψ2.*ψ2 .+ Dxψ2.*Dxψ2 .+ Dzψ2.*Dzψ2) +
            maximum(h1_out_x) +
            dρ*dz*sum(ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1 .+
                      Dρϕ1.*Dρϕ1 .+ Dρψv1.*Dρψv1 .+ Dρψ1.*Dρψ1 .+
                      Dzϕ1.*Dzϕ1 .+ Dzψv1.*Dzψv1 .+ Dzψ1.*Dzψ1) +
                      dt0*dz*sum(ψ1[1,:].*ψ1[1,:] .+ Dρψ1[1,:].*Dρψ1[1,:] .+
                                 Dzψ1[1,:].*Dzψ1[1,:])
    end

    tt, q, H1
    
end

# maximum times we double resolution
Nmax = 4

# base lowest resolution
Nρ = 17
Nz = 16

root_dir  = "/home/pmzag1/repos/model_CCE_CCM_public/examples/run_ccm/"
toy_model = "WH_WH_noise_t20_H1_amp/"

coarse_dir = joinpath(root_dir, toy_model, "data_$(Nρ)_$(Nz)")

x0 = h5read(coarse_dir * "/X.h5", "X")
ρ0 = h5read(coarse_dir * "/ρ.h5", "ρ")
dx0 = x0[2] - x0[1]
dρ0 = ρ0[2] - ρ0[1]
# timestep, CFL=0.25 (!VERIFY!)
dt0 = 0.25*minimum([dρ0, dx0])

for n in 0:1:Nmax

    @printf "n = %9d\n" n

    dir = joinpath(root_dir, toy_model, "data_$((Nρ-1)*2^n + 1)_$((Nz)*2^n)")
    
    tt, q, H1 = dev_t_func(dir, dt0)
    
    data_dir2 = joinpath(root_dir, toy_model, "dev_norms_exact")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "dev_norms_$(n).dat")
    open(outfile, "w") do io
        println(io, "# tt | q | H1 ")
        writedlm(io, [tt q H1])
    end

end
