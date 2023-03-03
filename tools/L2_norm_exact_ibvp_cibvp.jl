
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

    # full
    # for solution
    L2_cauchy_sol          = zeros(Nf)
    L2_char_sol            = zeros(Nf)
    # for given data
    L2_cauchy_given        = zeros(Nf)
    L2_char_given          = zeros(Nf)

    # temporary values for calculation of worldtube sums
    # out stands for outgoing vars, in for ingoing vars
    # characteristic
    l2_x_out     = zeros(length(ρ)) # sum for outgoing vars over all x
    l2_x_max_out = 0.0 # the max value of that sum, no necessarily xmax or xmin
    l2_xmin_out  = 0.0
    l2_xmax_in   = 0.0
    
    #cauchy
    l2_cauchy_ρmin_in  = 0.0
    l2_cauchy_ρmin_out = 0.0
    l2_cauchy_ρmax_in  = 0.0
    l2_cauchy_ρmax_out = 0.0

    # Initial data
    l2_char_u0   = 0.0
    l2_cauchy_t0 = 0.0
    
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

        if i==1
            println("computing initial data norm")
            l2_char_u0   = dx*dz*sum(ψ2.*ψ2)
            l2_cauchy_t0 = dρ*dz*sum(ϕ1.*ϕ1 + ψv1.*ψv1 + ψ1.*ψ1)
        end
        
        # sum the outgoing and ingoing norms to get the complete one
        # characteristic
        l2_x_out     += dt0*dz*sum(ϕ2.*ϕ2 + ψv2.*ψv2, dims=2) # sol
        l2_x_max_out  = maximum(l2_x_out) # sol
        l2_xmin_in   += dt0*dz*sum(ψ2[1,:].*ψ2[1,:]) # sol
        l2_xmin_out  += dt0*dz*sum(ϕ2[1,:].*ϕ2[1,:] + ψv2[1,:].*ψv2[1,:]) # given
        l2_xmax_in   += dt0*dz*sum(ψ2[end,:].*ψ2[end,:]) # given
        # cauchy
        l2_cauchy_ρmin_in  += dt0*dz*sum(ψ1[1,:].*ψ1[1,:]) # sol
        l2_cauchy_ρmax_out += dt0*dz*sum(ϕ1[end,:].*ϕ1[end,:] + ψv1[end,:].*ψv1[end,:]) # sol
        l2_cauchy_ρmax_in  += dt0*dz*sum(ψ1[end,:].*ψ1[end,:]) # given
        l2_cauchy_ρmin_out += dt0*dz*sum(ϕ1[1,:].*ϕ1[1,:] + ψv1[1,:].*ψv1[1,:]) # given
                        
        # full norms
        # solutions
        L2_cauchy_sol[it] = dρ*dz*sum(ϕ1.*ϕ1 + ψv1.*ψv1 + ψ1.*ψ1) +
            l2_cauchy_ρmin_in + l2_cauchy_ρmax_out
        L2_char_sol[it] = dx*dz*sum(ψ2.*ψ2) + l2_xmin_in + l2_x_max_out

        # given data
        L2_cauchy_given[it] = l2_cauchy_t0 + l2_cauchy_ρmin_out + l2_cauchy_ρmax_in
        L2_char_given[it] = l2_char_u0 + l2_xmin_out + l2_xmax_in

    end

    tt, L2_cauchy_sol, L2_char_sol, L2_cauchy_given, L2_char_given  
end

# maximum times we double resolution
Nmax = 4

# base lowest resolution
Nρ = 17
Nz = 16

root_dir  = "/home/thanasis/repos/model_CCE_CCM_public/examples/run_ibvp_cibvp/"
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
    
    tt, L2_cauchy_sol, L2_char_sol, L2_cauchy_given, L2_char_given = L2_t_func(dir, dt0)

    data_dir2 = joinpath(root_dir, toy_model, "norms_exact_cauchy_char")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "L2_$(n).dat")
    open(outfile, "w") do io
        println(io, "# tt | L2_cauchy_sol | L2_char_sol | L2_cauchy_given | L2_char_given")
        writedlm(io, [tt L2_cauchy_sol L2_char_sol L2_cauchy_given L2_char_given])
    end

end
