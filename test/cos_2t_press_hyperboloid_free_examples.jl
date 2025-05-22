# Pressurized hyperboloid entirely free.
# 
# Example introduced in
# @article{Lee2004,
#    author = {Lee, P. S. and Bathe, K. J.},
#    title = {Development of MITC isotropic triangular shell finite elements},
#    journal = {Computers & Structures},
#    volume = {82},
#    number = {11-12},
#    pages = {945-962},
#    ISSN = {0045-7949},
#    DOI = {10.1016/j.compstruc.2004.02.004},
#    year = {2004},
#    type = {Journal Article}
# }

module cos_2t_press_hyperboloid_free_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtools.MeshModificationModule: distortblock
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# Parameters:
E = 2.0e11 * phun("N")
nu = 1/3;
pressure = 1.0e6 * phun("Pa");
Length = 2.0 * phun("m");

# The hyperboloid axis is parallel to Y

function hyperbolic!(csmatout, XYZ, tangents, feid, qpid) 
    n = cross(tangents[:, 1], tangents[:, 2]) 
    n = n/norm(n)
    # r = vec(XYZ); r[2] = 0.0
    csmatout[:, 3] .= n
    csmatout[:, 2] .= (0.0, 1.0, 0.0)
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function computetrac!(forceout, XYZ, tangents, feid, qpid)
    r = vec(XYZ); r[2] = 0.0
    r .= vec(r)/norm(vec(r))
    theta = atan(r[3], r[1])
    n = cross(tangents[:, 1], tangents[:, 2]) 
    n = n/norm(n)
    forceout[1:3] = n*pressure*cos(2*theta)
    forceout[4:6] .= 0.0
    # @show dot(n, forceout[1:3])
    return forceout
end

function _execute(n = 8, thickness = Length/2/100, visualize = false, distortion = 0.0)
    formul = FinEtoolsFlexStructures.FEMMShellT3FFModule
    tolerance = Length/n/100
    fens, fes = distortblock(T3block, 90/360*2*pi, Length/2, n, n, distortion, distortion);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        R = sqrt(1 + y^2)
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0 0], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    # clamped edge perpendicular to Y
    # l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    # for i in [1,2,3,4,5,6]
    #     setebc!(dchi, l1, true, i)
    # end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    # nl = selectnode(fens; box = Float64[R R L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    
    fi = ForceIntensity(Float64, 6, computetrac!);
    F = distribloads(lfemm, geom0, dchi, fi, 2);
    
    # Solve
    solve_blocked!(dchi, K, F)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    strainenergy = 1/2 * U' * K * U
    @info "Strain Energy: $(round(strainenergy, digits = 9))"

    # Generate a graphical display of resultants
    ocsys = CSys(3, 3, hyperbolic!)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    visualize && vtkwrite("cos_2t_press_hyperboloid_free-$(n)-$(thickness)-$(distortion)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("n$nc", fld.values))
    end
    visualize && vtkwrite("cos_2t_press_hyperboloid_free-$(n)-$(thickness)-$(distortion)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end
    visualize && vtkwrite("cos_2t_press_hyperboloid_free-$(n)-$(thickness)-$(distortion)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])

    # Visualization
    if visualize
        scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end

    return strainenergy
end

function test_convergence(thicknessmult = 1/100, distortion = 0.0)
    println("--------------------------------------------------------------------------------")
    @info "Pressurized hyperboloid of rotation, free ends, thickness multiplier = $(thicknessmult), distortion = $(distortion)"
    results = []
    ns = [16, 32, 64,]
    for n in ns
        push!(results, _execute(n, Length/2*thicknessmult, false, 2*distortion/n))
    end
    return ns, results
end

end # module

using .cos_2t_press_hyperboloid_free_examples


let

    for distortion in [2.0, 0.0]
        ns, results100 = cos_2t_press_hyperboloid_free_examples.test_convergence(1/100, distortion)
        ns, results1000 = cos_2t_press_hyperboloid_free_examples.test_convergence(1/1000, distortion)
        ns, results10000 = cos_2t_press_hyperboloid_free_examples.test_convergence(1/10000, distortion)
        ns, results100000 = cos_2t_press_hyperboloid_free_examples.test_convergence(1/100000, distortion)

        @show q1, q2, q3 = results100
        @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
        @show q1, q2, q3 = results1000
        @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
        @show q1, q2, q3 = results10000
        @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
        @show q1, q2, q3 = results100000
        @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

    end
end