
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


function L2_cmh_t(dir_c, dir_m, dir_h)
    # ρ grid
    ρc = h5read(dir_c * "/ρ.h5", "ρ")
    ρm = h5read(dir_m * "/ρ.h5", "ρ")
    ρh = h5read(dir_h * "/ρ.h5", "ρ")
    # x grid
    xc = h5read(dir_c * "/X.h5", "X")
    xm = h5read(dir_m * "/X.h5", "X")
    xh = h5read(dir_h * "/X.h5", "X")
    # z grid
    zc = h5read(dir_c * "/z.h5", "z")
    zm = h5read(dir_m * "/z.h5", "z")
    zh = h5read(dir_h * "/z.h5", "z")

    # make sure that we can inject points from the medium and high resolution
    # grids in the coarse grid without interpolation
    @assert ρc ≈ ρm[1:2:end] ≈ ρh[1:4:end]
    @assert xc ≈ xm[1:2:end] ≈ xh[1:4:end]
    @assert zc ≈ zm[1:2:end] ≈ zh[1:4:end]

    # dρ of coarse resolution
    dρc = ρc[2] - ρc[1] 
    # dx of coarse resolution
    dxc = xc[2] - xc[1] 
    # dz of coarse resolution
    dzc = zc[2] - zc[1]

    # list all available iterations (and corresponding files)
    (its_c, all_filenames_c) = list_h5_files(dir_c, prefix="data_")
    (its_m, all_filenames_m) = list_h5_files(dir_m, prefix="data_")
    (its_h, all_filenames_h) = list_h5_files(dir_h, prefix="data_")

    # we want to compare the common timesteps. we assume here that there is a
    # factor of 2 between the lowest resolution and the corresponding higher
    # resolution ones
    filenames_c = all_filenames_c[:]
    filenames_m = all_filenames_m[:]# if all timesteps are saved[1:2:end]
    filenames_h = all_filenames_h[:]#[1:4:end]

    Nf = length(filenames_c)
    @assert length(filenames_m) == length(filenames_h) == Nf

    tt     = zeros(Nf)
    L2_cmt = zeros(Nf)
    L2_mht = zeros(Nf)
    
    # initiate the grid function for the outgoing norm
    l2_out_x_cmt = zeros(length(xc))
    l2_out_x_mht = zeros(length(xc))

    l2_cauchy_ρmin_in_cmt = 0.0
    l2_cauchy_ρmin_in_mht = 0.0
    
    for it in ProgressBar(1:Nf)
        
        file_c = filenames_c[it]
        file_m = filenames_m[it]
        file_h = filenames_h[it]

        ψ1c  = h5read(file_c, "ψ1")
        ψv1c = h5read(file_c, "ψv1")
        ϕ1c  = h5read(file_c, "ϕ1")
        ψ2c  = h5read(file_c, "ψ2")
        ψv2c = h5read(file_c, "ψv2")
        ϕ2c  = h5read(file_c, "ϕ2")
        tc  = h5readattr(file_c, "/")["time"]

        ψ1m  = h5read(file_m, "ψ1")
        ψv1m = h5read(file_m, "ψv1")
        ϕ1m  = h5read(file_m, "ϕ1")
        ψ2m  = h5read(file_m, "ψ2")
        ψv2m = h5read(file_m, "ψv2")
        ϕ2m  = h5read(file_m, "ϕ2")
        tm  = h5readattr(file_m, "/")["time"]

        ψ1h  = h5read(file_h, "ψ1")
        ψv1h = h5read(file_h, "ψv1")
        ϕ1h  = h5read(file_h, "ϕ1")
        ψ2h  = h5read(file_h, "ψ2")
        ψv2h = h5read(file_h, "ψv2")
        ϕ2h  = h5read(file_h, "ϕ2")
        th  = h5readattr(file_h, "/")["time"]

        # make sure we're comparing the same timestep
        @assert tc ≈ tm ≈ th

        # compute the differences between resolutions (projected onto the coarsest grid)
        ψ1cm  = ψ1c  .-  ψ1m[1:2:end, 1:2:end]
        ψ2cm  = ψ2c  .-  ψ2m[1:2:end, 1:2:end]
        ψv1cm = ψv1c .- ψv1m[1:2:end, 1:2:end]
        ψv2cm = ψv2c .- ψv2m[1:2:end, 1:2:end]
        ϕ1cm  = ϕ1c  .-  ϕ1m[1:2:end, 1:2:end]
        ϕ2cm  = ϕ2c  .-  ϕ2m[1:2:end, 1:2:end]
        ψ1mh  =  ψ1m[1:2:end, 1:2:end]  .- ψ1h[1:4:end, 1:4:end]
        ψ2mh  =  ψ2m[1:2:end, 1:2:end]  .- ψ2h[1:4:end, 1:4:end]
        ψv1mh = ψv1m[1:2:end, 1:2:end] .- ψv1h[1:4:end, 1:4:end]
        ψv2mh = ψv2m[1:2:end, 1:2:end] .- ψv2h[1:4:end, 1:4:end]
        ϕ1mh  =  ϕ1m[1:2:end, 1:2:end]  .- ϕ1h[1:4:end, 1:4:end]
        ϕ2mh  =  ϕ2m[1:2:end, 1:2:end]  .- ϕ2h[1:4:end, 1:4:end]

        tt[it]     = tc

        # define the timestep; needed for the sum in u
        # give the same cfl condition as in the runs
        dt0 = 0.25*minimum([dρc, dxc])

        # get the characteristi outgoing part of the norm
        l2_out_x_cmt += dt0*dzc*sum( ψv2cm.*ψv2cm .+ ϕ2cm.*ϕ2cm, dims=2)
        l2_out_x_mht += dt0*dzc*sum( ψv2mh.*ψv2mh .+ ϕ2mh.*ϕ2mh, dims=2)
        # cauchy in at left boundary
        l2_cauchy_ρmin_in_cmt += dt0*dzc*sum(ψ1cm[1,:].*ψ1cm[1,:])
        l2_cauchy_ρmin_in_mht += dt0*dzc*sum(ψ1mh[1,:].*ψ1mh[1,:])
        
        # sum the outgoing and ingoing norms to get the complete one
        L2_cmt[it] = dxc*dzc*sum(ψ2cm.*ψ2cm) + maximum(l2_out_x_cmt) +
            dρc*dzc*sum(ϕ1cm.*ϕ1cm .+ ψv1cm.*ψv1cm .+ ψ1cm.*ψ1cm) +
            l2_cauchy_ρmin_in_cmt
        
        L2_mht[it] = dxc*dzc*sum(ψ2mh.*ψ2mh) + maximum(l2_out_x_mht) +
            dρc*dzc*sum(ϕ1mh.*ϕ1mh .+ ψv1mh.*ψv1mh .+ ψ1mh.*ψ1mh) +
            l2_cauchy_ρmin_in_mht
        
    end
    
    tt, L2_cmt, L2_mht
end


# maximum times we double resolution
Nmax = 3

# base lowest resolution
Nx = 17
Nz = 16

root_dir  = "/home/thanasis/repos/model_CCE_CCM_public/examples/run_ccm/"
toy_model = "SYMH_B1_WH_B2_smooth_t20/"
#"advect_noise_t12_0.25/"

# we need 3 different resolutions to build the L2 norm that is used in the self
# convergence ratio
for n in 0:1:Nmax-2
    @printf "n = %9d\n" n

    dir_c = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^n + 1)_$(Nz*2^n)")
    dir_m = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))")
    dir_h = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))")

    tt, L2_cmt, L2_mht = L2_cmh_t(dir_c, dir_m, dir_h)

    data_dir2 = joinpath(root_dir, toy_model, "norms_self")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "L2_$(n)_$(n+1)_c$(n).dat")
    open(outfile, "w") do io
        println(io, "# t      |      L2")
        writedlm(io, [tt L2_cmt])
    end

    outfile  = joinpath(data_dir2, "L2_$(n+1)_$(n+2)_c$(n).dat")
    open(outfile, "w") do io
        println(io, "# t      |      L2")
        writedlm(io, [tt L2_mht])
    end
end
