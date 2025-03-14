using VortexLattice
using Plots


root_chord = 2.2
span = 15
ns = 15

function define_elliptical(root_chord, span, ns)
    chord_distribution(y) = root_chord * sqrt(1-(y/(span/2))^2)

    yle = range(0, stop=span/2, length=ns)
    xle = zeros(ns)
    zle = zeros(ns)
    chord = [chord_distribution(y) for y in yle]
    
    return yle, xle, zle, chord
end


function elliptical(root_chord, span, ns)
    yle, xle, zle, chord = define_elliptical(root_chord, span, ns)
    theta = fill(2.0*pi/180, ns)
    phi = fill(0.0, ns)
    fc = fill((xc) -> 0, ns) #camberline function for each section

    
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # declare symmetry
    symmetric = true

    # construct surface
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

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
    write_vtk("elliptical-wing", system)
end

elliptical(root_chord, span, ns)