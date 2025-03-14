using VortexLattice
using Plots


root_chord = 2.2
span = 15.0

function define_elliptical(root_chord, span, ns)
    chord_distribution(y) = root_chord * sqrt(1-(y/(span/2))^2) #defines chord distribution by modification of ellipse equation

    yle = range(0, stop=span/2, length=ns) #defines y component of wing as a range from 0 to half the span, divided into ns sections
    xle = [root_chord / 2 * (1-sqrt(1-(2*y/span)^2)) for y in yle] # rearranged elliptical equation to define leading edge position in x component
    zle = zeros(ns) #defines z component of wing as a vector of zeros, same length as ns
    chord = [chord_distribution(y) for y in yle]   #defines chord length at each section by applying chord_distribution function to each element in yle vector
    
    return yle, xle, zle, chord #returns these values to be used in VLM calculations
end


function elliptical(root_chord, span, ns)
    yle, xle, zle, chord = define_elliptical(root_chord, span, ns)
    theta = fill(2.0*pi/180, ns)
    phi = fill(0.0, ns)
    fc = fill((xc) -> 0, ns) #camberline function for each section

    
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

    Sref = pi/4 * root_chord * span #reference area for the wing, calculated as a quarter of the product of root chord and span
    cref = 2.2
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
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
    efficiency = CL / CD

    #write_vtk("elliptical-wing", system)

    return efficiency
end

ns_values = 2:1:50
efficiencies = [elliptical(root_chord, span, ns) for ns in ns_values]


plot(ns_values, efficiencies, xlabel="Number of sections (ns)", ylabel="Efficiency", legend =false)