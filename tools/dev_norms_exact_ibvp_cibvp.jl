
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
function norms_t_func(dir, dt0)
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
    # for solution
    q_cauchy_sol     = zeros(Nf)
    q_char_sol       = zeros(Nf)
    H1_cauchy_sol    = zeros(Nf)
    H1_char_sol      = zeros(Nf)
    # for given data
    q_cauchy_given   = zeros(Nf)
    q_char_given     = zeros(Nf)
    H1_cauchy_given  = zeros(Nf)
    H1_char_given    = zeros(Nf)

    # temporary values for calculation of worldtube sums
    # out stands for outgoing vars, in for ingoing vars
    # characteristic
    q_x_out       = zeros(length(ρ)) # sum for outgoing vars over all x
    q_x_max_out   = 0.0 # the max value of that sum, no necessarily xmax or xmin
    h1_x_out      = zeros(length(ρ)) # sum for outgoing vars over all x
    h1_x_max_out  = 0.0 # the max value of that sum
    q_xmin_out    = 0.0
    q_xmin_in     = 0.0
    q_xmax_in     = 0.0
    h1_xmin_in    = 0.0
    h1_xmin_out   = 0.0
    h1_xmax_out   = 0.0
    h1_xmax_in    = 0.0
    
    #cauchy
    q_cauchy_ρmin_in   = 0.0
    q_cauchy_ρmin_out  = 0.0
    q_cauchy_ρmax_in   = 0.0
    q_cauchy_ρmax_out  = 0.0
    h1_cauchy_ρmin_in  = 0.0
    h1_cauchy_ρmin_out = 0.0
    h1_cauchy_ρmax_in  = 0.0
    h1_cauchy_ρmax_out = 0.0

    # Initial data
    q_char_u0   = 0.0
    q_cauchy_t0 = 0.0
    h1_char_u0   = 0.0
    h1_cauchy_t0 = 0.0
    
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
        
        # time
        tt[it]  = t        

        if it==1
            q_char_u0    = dx*dz*sum(ψ2.*ψ2)
            q_cauchy_t0  = dρ*dz*sum(ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1 .+ Dzϕ1.*Dzϕ1)
            h1_char_u0   = dx*dz*sum(ψ2.*ψ2 .+ Dxψ2.*Dxψ2 .+ Dzψ2.*Dzψ2)
            h1_cauchy_t0 = dρ*dz*sum(ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1 .+
                                     Dρϕ1.*Dρϕ1 .+ Dρψv1.*Dρψv1 .+ Dρψ1.*Dρψ1 .+
                                     Dzϕ1.*Dzϕ1 .+ Dzψv1.*Dzψv1 .+ Dzψ1.*Dzψ1)
        end

        # sum the outgoing and ingoing norms to get the complete one
        # characteristic
        q_x_out     += dt0*dz*sum(ϕ2.*ϕ2 + ψv2.*ψv2 + Dzϕ2.*Dzϕ2, dims=2) # sol
        q_x_max_out  = maximum(q_x_out) # sol
        q_xmin_in   += dt0*dz*sum(ψ2[1,:].*ψ2[1,:]) # sol
        q_xmin_out  += dt0*dz*sum(ϕ2[1,:].*ϕ2[1,:] + ψv2[1,:].*ψv2[1,:] +
                                  Dzϕ2[1,:].*Dzϕ2[1,:]) # given
        q_xmax_in   += dt0*dz*sum(ψ2[end,:].*ψ2[end,:]) # given
        # cauchy
        q_cauchy_ρmin_in  += dt0*dz*sum(ψ1[1,:].*ψ1[1,:]) # sol
        q_cauchy_ρmax_out += dt0*dz*sum(ϕ1[end,:].*ϕ1[end,:] + ψv1[end,:].*ψv1[end,:] +
                                        Dzϕ1[end,:].*Dzϕ1[end,:]) # sol
        q_cauchy_ρmax_in  += dt0*dz*sum(ψ1[end,:].*ψ1[end,:]) # given
        q_cauchy_ρmin_out += dt0*dz*sum(ϕ1[1,:].*ϕ1[1,:] + ψv1[1,:].*ψv1[1,:] +
                                        Dzϕ1[1,:].*Dzϕ1[1,:]) # given
        # characteristic
        h1_x_out     += dt0*dz*sum(ϕ2.*ϕ2 + ψv2.*ψv2 +
                                   Dxϕ2.*Dxϕ2 + Dxψv2.*Dxψv2 +
                                   Dzϕ2.*Dzϕ2 + Dzψv2.*Dzψv2, dims=2) # sol
        h1_x_max_out  = maximum(h1_x_out) # sol
        h1_xmin_in   += dt0*dz*sum(ψ2[1,:].*ψ2[1,:] +
                                   Dxψ2[1,:].*Dxψ2[1,:] +
                                   Dzψ2[1,:].*Dzψ2[1,:]) # sol
        h1_xmin_out  += dt0*dz*sum(ϕ2[1,:].*ϕ2[1,:] + ψv2[1,:].*ψv2[1,:] +
                                   #Dxϕ2[1,:].*Dxϕ2[1,:] + Dxψv2[1,:].*Dxψv2[1,:]+ # involves part of sol to find it; not controlled
                                   Dzϕ2[1,:].*Dzϕ2[1,:] + Dzψv2[1,:].*Dzψv2[1,:]) # given
        h1_xmax_in   += dt0*dz*sum(ψ2[end,:].*ψ2[end,:] +
                                   #Dxψ2[end,:].*Dxψ2[end,:] + # involves part of sol to find it; not controlled
                                   Dzψ2[end,:].*Dzψ2[end,:]) # given
        # cauchy
        h1_cauchy_ρmin_in  += dt0*dz*sum(ψ1[1,:].*ψ1[1,:] +
                                         Dρψ1[1,:].*Dρψ1[1,:] +
                                         Dzψ1[1,:].*Dzψ1[1,:]) # sol
        h1_cauchy_ρmax_out += dt0*dz*sum(
            ϕ1[end,:].*ϕ1[end,:] + ψv1[end,:].*ψv1[end,:] +
            Dρϕ1[end,:].*Dρϕ1[end,:] + Dρψv1[end,:].*Dρψv1[end,:] +
            Dzϕ1[end,:].*Dzϕ1[end,:] + Dzψv1[end,:].*Dzψv1[end,:]) # sol
        h1_cauchy_ρmax_in  += dt0*dz*sum(ψ1[end,:].*ψ1[end,:] +
                                         #Dρψ1[end,:].*Dρψ1[end,:] + # involves part of sol to find it; not controlled
                                         Dzψ1[end,:].*Dzψ1[end,:]) # given
        h1_cauchy_ρmin_out += dt0*dz*sum(
            ϕ1[1,:].*ϕ1[1,:] + ψv1[1,:].*ψv1[1,:] +
            #Dρϕ1[1,:].*Dρϕ1[1,:] + Dρψv1[1,:].*Dρψv1[1,:] + # involves part of sol to find it; not controlled
            Dzϕ1[1,:].*Dzϕ1[1,:] + Dzψv1[1,:].*Dzψv1[1,:]) # given

        # full norms
        q_cauchy_sol[it] = dρ*dz*sum(ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1 .+ Dzϕ1.*Dzϕ1) +
            q_cauchy_ρmin_in + q_cauchy_ρmax_out

        q_char_sol[it] = dx*dz*sum(ψ2.*ψ2) + q_xmin_in + q_x_max_out

        H1_cauchy_sol[it] = dρ*dz*sum(
            ϕ1.*ϕ1 .+ ψv1.*ψv1 .+ ψ1.*ψ1 .+
            Dρϕ1.*Dρϕ1 .+ Dρψv1.*Dρψv1 .+ Dρψ1.*Dρψ1 .+
            Dzϕ1.*Dzϕ1 .+ Dzψv1.*Dzψv1 .+ Dzψ1.*Dzψ1) +
            h1_cauchy_ρmin_in + h1_cauchy_ρmax_out

        H1_char_sol[it] = dx*dz*sum(ψ2.*ψ2 .+ Dxψ2.*Dxψ2 .+ Dzψ2.*Dzψ2) +
            h1_xmin_in + h1_x_max_out

        # given data norms
        q_cauchy_given[it]  = q_cauchy_t0 + q_cauchy_ρmin_out + q_cauchy_ρmax_in
        q_char_given[it]    = q_char_u0 + q_xmin_out + q_xmax_in
        H1_cauchy_given[it] = h1_cauchy_t0 + h1_cauchy_ρmin_out + h1_cauchy_ρmax_in
        H1_char_given[it]   = h1_char_u0 + h1_xmin_out + h1_xmax_in
            
    end

    tt, q_cauchy_sol, q_char_sol, H1_cauchy_sol, H1_char_sol, q_cauchy_given, q_char_given, H1_cauchy_given, H1_char_given

end

# maximum times we double resolution
Nmax = 4

# base lowest resolution
Nρ = 17
Nz = 16

root_dir  = "/home/pmzag1/repos/model_CCE_CCM_public/examples/run_cce/"
toy_model = "SYMH_WH_noise_t20_q_amp/"

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
    
    
    tt, q_cauchy_sol, q_char_sol, H1_cauchy_sol, H1_char_sol, q_cauchy_given, q_char_given, H1_cauchy_given, H1_char_given = norms_t_func(dir, dt0)

    data_dir2 = joinpath(root_dir, toy_model, "dev_norms_exact_cauchy_char")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "norms_$(n).dat")
    open(outfile, "w") do io
        println(io, "# tt | q_cauchy_sol | q_char_sol | H1_cauchy_sol | H1_char_sol | q_cauchy_given | q_char_given | H1_cauchy_given | H1_char_given")
        writedlm(io, [tt q_cauchy_sol q_char_sol H1_cauchy_sol H1_char_sol q_cauchy_given q_char_given H1_cauchy_given H1_char_given])
    end

end
