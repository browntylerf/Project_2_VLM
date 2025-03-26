using VortexLattice
using Plots

function varying_alpha()
    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [0.0*pi/180, 0.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) #camberline function for each section

    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha_range = -20:0.5:20
    alphas = alpha_range .*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]

    CL_values = []

    for alpha in alphas
        fs = Freestream(Vinf, alpha, beta, Omega)

        # construct surface
        grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=mirror)

        # create vector containing all grids and ratios
        grids = [grid]
        ratios = [ratio]

     # create the system object
        system = System(grids; ratios)

        # perform steady state analysis
        steady_analysis!(system, ref, fs; symmetric=symmetric)

        # retrieve near-field forces
        CF, CM = body_forces(system; frame=Wind())

        # perform far-field analysis
        CDiff = far_field_drag(system)

        CD, CY, CL = CF
        Cl, Cm, Cn = CM
        push!(CL_values, CL)
        write_vtk("varying_alpha", system)
    end
    plot(alpha_range, CL_values, xlabel="Alpha (degrees)", ylabel="CL", legend=false)
    savefig("Varying_Alpha.png")
end

varying_alpha()