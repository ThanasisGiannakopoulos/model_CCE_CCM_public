
using HDF5
using DelimitedFiles
using Printf
using ProgressBars

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
function L2_t_func(dir, dt0)
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

    #L2_t         = zeros(Nf)
    # full
    L2_cauchy        = zeros(Nf)
    L2_char          = zeros(Nf)
    #characteristic
    l2_out_x     = zeros(length(ρ))
    l2_out_x_max = zeros(Nf)
    
    # temporary values for calculation
    # characteristic
    l2_out_end_temp = 0.0
    #cauchy
    l2_cauchy_in_1_temp = 0.0
    l2_cauchy_in_ρ_temp = zeros(length(ρ))
    
    for it in ProgressBar(1:Nf)
        file = filenames[it]

        # cauchy
        ϕ1  = h5read(file, "ϕ1")
        ψv1 = h5read(file, "ψv1")
        ψ1  = h5read(file, "ψ1")
        # charteristic
        ϕ2  = h5read(file, "ϕ2")
        ψv2 = h5read(file, "ψv2")
        ψ2  = h5read(file, "ψ2")
        # time
        t  = h5readattr(file, "/")["time"]
        
        # time
        tt[it]  = t

        # sum the outgoing and ingoing norms to get the complete one
        # characteristic
        l2_out_x           += dt0*dz*sum(ϕ2.*ϕ2 + ψv2.*ψv2, dims=2)
        l2_out_x_max[it]    = maximum(l2_out_x)
        
        # full norms
        L2_cauchy[it] = dρ*dz*sum(ϕ1.*ϕ1 + ψv1.*ψv1 + ψ1.*ψ1) +
            dt0*dz*sum(ψ1[1,:].*ψ1[1,:]) +
            dt0*dz*sum(ϕ1[end,:].*ϕ1[end,:] .+ ψv1[end,:].*ψv1[end,:])

        L2_char[it]   = dx*dz*sum(ψ2.*ψ2) +
            dt0*dz*sum(ψ2[1,:].*ψ2[1,:]) + l2_out_x_max[it]
       
    end

    tt, L2_cauchy, L2_char 
end

# maximum times we double resolution
Nmax = 4

# base lowest resolution
Nρ = 17
Nz = 16

root_dir  = "/home/thanasis/repos/model_CCE_CCM_public/examples/run_cce/"
toy_model = "SYMH_SYMH_noise_t20_L2_amp/"

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
    
    tt, L2_cauchy, L2_char = L2_t_func(dir, dt0)

    data_dir2 = joinpath(root_dir, toy_model, "norms_exact_cauchy_char")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "L2_$(n).dat")
    open(outfile, "w") do io
        println(io, "# tt | L2_cauchy | L2_char")
        writedlm(io, [tt L2_cauchy L2_char])
    end

end
