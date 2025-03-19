using VortexLattice
using Plots

function volume_ratios()
    # wing
    xle = [0.0, 0.0]
    yle = [0.0, 5.0]
    zle = [0.0, 0.0]
    chord = [1.0, 1.0]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true

    Sref = 5.0 
    cref = 1.0
    bref = 5.0
    rref = [0.5, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 0.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    span_h_range = range(1.0, 3.0, 10)
    span_v_range = range(1.0, 3.0, 10)

    horiz_coeff = Float64[]
    vert_coeff = Float64[]
    Cla_horiz = Float64[]
    Cma_horiz = Float64[]
    Cna_horiz = Float64[]
    Cla_vert = Float64[]
    Cma_vert = Float64[]
    Cna_vert = Float64[]

    for span_h in span_h_range
        # horizontal stabilizer
        xle_h = [0.0, 0.0]
        yle_h = [0.0, span_h]
        zle_h = [0.0, 0.0]
        chord_h = [0.7, 0.7]
        theta_h = [0.0, 0.0]
        phi_h = [0.0, 0.0]
        fc_h = fill((xc) -> 0, 2) # camberline function for each section
        ns_h = 6
        nc_h = 3
        spacing_s_h = Uniform()
        spacing_c_h = Uniform()
        mirror_h = true


        # vertical stabilizer
        xle_v = [0.0, 0.0]
        yle_v = [0.0, 0.0]
        zle_v = [0.0, 1.0]
        chord_v = [0.7, 0.7]
        theta_v = [0.0, 0.0]
        phi_v = [0.0, 0.0]
        fc_v = fill((xc) -> 0, 2) # camberline function for each section
        ns_v = 5
        nc_v = 3
        spacing_s_v = Uniform()
        spacing_c_v = Uniform()
        mirror_v = false

        symmetric = [false, false, false]

        # generate surface panels for wing
        wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
            mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)
    
        # generate surface panels for horizontal tail
        hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
            mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
        VortexLattice.translate!(hgrid, [4.0, 0.0, 0.0])
    
        # generate surface panels for vertical tail
        vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
            mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
        VortexLattice.translate!(vgrid, [4.0, 0.0, 0.0])
    
        # now set normal vectors manually
        ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 0.0*pi/180)
    
        grids = [wgrid, hgrid, vgrid]
        ratios = [wratio, hratio, vratio]
        surface_id = [1, 2, 3]
    
        system = System(grids; ratios)
    
        # overwrite normal vector for each wing panel
        for i = 1:length(system.surfaces[1])
            system.surfaces[1][i] = set_normal(system.surfaces[1][i], ncp)
        end
    
        steady_analysis!(system, ref, fs; symmetric=symmetric, surface_id=surface_id, derivatives=true)
    
        CF, CM = body_forces(system; frame=Stability())
    
        CDiff = far_field_drag(system)
    
        CD, CY, CL = CF
        Cl, Cm, Cn = CM
    
        dCF, dCM = stability_derivatives(system)
        
        Cla, Cma, Cna = dCM[:alpha]
        Clb, Cmb, Cnb = dCM[:beta]

        write_vtk("wing-tail", system)

        tail_volume_coeff_horiz = (span_h*chord_h[1]*4)/(Sref * cref) #4 is the distance between the aerodynamic center of the wing and the horizontal tail

        push!(horiz_coeff, tail_volume_coeff_horiz)
        push!(Cla_horiz, Cma) #longitudinal stability
        push!(Cma_horiz, Cnb) #yaw stability
        push!(Cna_horiz, Clb) #roll stability
    end


    for span_v in span_v_range
        # horizontal stabilizer
        xle_h = [0.0, 0.0]
        yle_h = [0.0, 1.0]
        zle_h = [0.0, 0.0]
        chord_h = [0.7, 0.7]
        theta_h = [0.0, 0.0]
        phi_h = [0.0, 0.0]
        fc_h = fill((xc) -> 0, 2) # camberline function for each section
        ns_h = 6
        nc_h = 3
        spacing_s_h = Uniform()
        spacing_c_h = Uniform()
        mirror_h = true


        # vertical stabilizer
        xle_v = [0.0, 0.0]
        yle_v = [0.0, 0.0]
        zle_v = [0.0, span_v]
        chord_v = [0.7, 0.7]
        theta_v = [0.0, 0.0]
        phi_v = [0.0, 0.0]
        fc_v = fill((xc) -> 0, 2) # camberline function for each section
        ns_v = 5
        nc_v = 3
        spacing_s_v = Uniform()
        spacing_c_v = Uniform()
        mirror_v = false

        symmetric = [false, false, false]

        # generate surface panels for wing
        wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
            mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)
    
        # generate surface panels for horizontal tail
        hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
            mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
        VortexLattice.translate!(hgrid, [4.0, 0.0, 0.0])
    
        # generate surface panels for vertical tail
        vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
            mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
        VortexLattice.translate!(vgrid, [4.0, 0.0, 0.0])
    
        # now set normal vectors manually
        ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 0.0*pi/180)
    
        grids = [wgrid, hgrid, vgrid]
        ratios = [wratio, hratio, vratio]
        surface_id = [1, 2, 3]
    
        system = System(grids; ratios)
    
        # overwrite normal vector for each wing panel
        for i = 1:length(system.surfaces[1])
            system.surfaces[1][i] = set_normal(system.surfaces[1][i], ncp)
        end
    
        steady_analysis!(system, ref, fs; symmetric=symmetric, surface_id=surface_id, derivatives=true)
    
        CF, CM = body_forces(system; frame=Stability())
    
        CDiff = far_field_drag(system)
    
        CD, CY, CL = CF
        Cl, Cm, Cn = CM
    
        dCF, dCM = stability_derivatives(system)
        
        Cla, Cma, Cna = dCM[:alpha]
        Clb, Cmb, Cnb = dCM[:beta]
        

        write_vtk("wing-tail", system)

        tail_volume_coeff_vert = (span_v*chord_v[1]*4)/(Sref * cref) #4 is the distance between the aerodynamic center of the wing and the vertical tail

        push!(vert_coeff, tail_volume_coeff_vert)
        push!(Cla_vert, Cma) #longitudinal stability
        push!(Cma_vert, Cnb) #yaw stability
        push!(Cna_vert, Clb) #roll stability
    end
    print(Cma_vert)
end

function avl_normal_vector(vector, angle)
    normal_vector = [vector[1] * cos(angle), vector[2] * cos(angle), vector[3] * sin(angle)]
    return normal_vector
end

volume_ratios()

