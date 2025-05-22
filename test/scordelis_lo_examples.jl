"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

Problem description

The physical basis of the problem is a deeply arched roof supported only
by diaphragms at its curved edges (an aircraft hangar), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight (downwardly directed)
load.  
"""
module scordelis_lo_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

# The "reference" solution for the vertical deflection at the midpoint of the
# free edge. Note that this number is most likely wrong, or at least inaccurate (for shear flexible shells).
analyt_sol = -0.3024 * phun("ft")

# Parameters:
E = 4.32e8 * phun("lbf/ft^2")
nu = 0.0
thickness = 0.25 * phun("ft") 
R = 25.0 * phun("ft")
L = 50.0 * phun("ft")
q = 90.0 * phun("lbf/ft^2")

cylindrical!(csmatout, XYZ, tangents, fe_label, qpid) = begin
    r = vec(XYZ)
    r[2] = 0.0
    r[3] += R
    csmatout[:, 3] .= vec(r) / norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _execute(n = 8, visualize = true)
    formul = FEMMShellT3FFModule

    tolerance = R / n / 10
    # fens, fes = T3blockrand(40/360*2*pi,L/2,n,n);
    fens, fes = T3block(40 / 360 * 2 * pi, L / 2, n, n, :a)
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        y = fens.xyz[i, 2]
        fens.xyz[i, :] .= (R * sin(a), y, R * (cos(a) - 1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1, 3, 5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L / 2 L / 2 -Inf Inf], inflate = tolerance)
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi)

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi)

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40 / 360 * 2 * pi) * R sin(40 / 360 * 2 * pi) * R L / 2 L / 2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(Float64[0, 0, -q, 0, 0, 0])
    F = distribloads(lfemm, geom0, dchi, fi, 3)

    # Solve
    solve_blocked!(dchi, K, F)
    result = dchi.values[nl, 3][1]
    @info "Deflection: $(result / phun("ft")) [ft]"

    # Visualization
    if visualize
        # Generate a graphical display of resultants
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("m$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            push!(scalars, ("em$nc", fld.values))
        end
        vtkwrite("scordelis_lo_examples-$(n)-m.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:3
            fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("n$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            push!(scalars, ("en$nc", fld.values))
        end
        vtkwrite("scordelis_lo_examples-$(n)-n.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        scalars = []
        for nc in 1:2
            fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("q$nc", fld.values))
            fld = elemfieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            push!(scalars, ("eq$nc", fld.values))
        end
        vtkwrite("scordelis_lo_examples-$(n)-q.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
        # Generate a graphical display of the displacements and rotations
        vtkwrite("scordelis_lo_examples-$(n)-uur.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3]), ("ur", dchi.values[:, 4:6])])
    end

    result
end

function test_convergence(ns = [4, 6, 8, 16, 32, 64, 128, 256], visualize = false)
    println("--------------------------------------------------------------------------------")
    @info "The barrel vault (Scordelis-Lo) roof"
    results = []
    for n in ns
        v = _execute(n, visualize)
        push!(results, v)
    end
    return ns, results
end

end # module

using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol
using .scordelis_lo_examples


ns, results = scordelis_lo_examples.test_convergence()

@show ns
@show results ./ phun("ft")
q1, q2, q3 = results[end-2:end] ./ phun("ft")
@show     qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

nothing