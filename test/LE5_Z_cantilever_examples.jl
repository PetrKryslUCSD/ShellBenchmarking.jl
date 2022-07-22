"""
MODEL DESCRIPTION

Z-section cantilever under torsional loading.

Linear elastic analysis, Young's modulus = 210 GPa, Poisson's ratio = 0.3.

All displacements are fixed at X=0.

Torque of 1.2 MN-m applied at X=10. The torque is applied by two 
uniformly distributed shear loads of 0.6 MN at each flange surface.

Objective of the analysis is to compute the axial stress at X = 2.5 from fixed end.

NAFEMS REFERENCE SOLUTION

Axial stress at X = 2.5 from fixed end (point A) at the midsurface is -108 MPa.

Note: the NAFEMS solution is likely to be inaccurate
"""
module LE5_Z_cantilever_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT35
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function zcant!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)  
    r = vec(XYZ); 
    cross3!(r, view(tangents, :, 1), view(tangents, :, 2))
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 1] .= (1.0, 0.0, 0.0)
    cross3!(view(csmatout, :, 2), view(csmatout, :, 3), view(csmatout, :, 1))
    return csmatout
end


function _execute(input = "nle5xf3c.inp", nrefs = 0, visualize = true)
    formul = FEMMShellT3FFModule
    E = 210 * phun("GPa");
    nu = 0.3;
    L = 10.0 * phun("m");
    d = 1.0 * phun("m")
    thickness = 0.1 * phun("m")
    S = 0.6 * phun("mega*N")

    tolerance = thickness/1000
    output = import_ABAQUS(joinpath(dirname(@__FILE__()), input))
    fens = output["fens"]
    fes = output["fesets"][1]

    connected = findunconnnodes(fens, fes);
    fens, new_numbering = compactnodes(fens, connected);
    fes = renumberconn!(fes, new_numbering);

    for r in 1:nrefs
        fens, fes = T3refine(fens, fes)
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, zcant!)
    
    # Report
    @info "Mesh: $input, nrefs = $nrefs"

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
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);
    @show dchi.nfreedofs

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);
    
    # Load
    # nl = selectnode(fens; box = Float64[L L d d 0 0], inflate = tolerance)
    # loadbdry1 = FESetP1(reshape(nl, 1, 1))
    # lfemm1 = FEMMBase(IntegDomain(loadbdry1, PointRule()))
    # fi1 = ForceIntensity(FFlt[0, 0, +S, 0, 0, 0]);
    # nl = selectnode(fens; box = Float64[L L -d -d 0 0], inflate = tolerance)
    # loadbdry2 = FESetP1(reshape(nl, 1, 1))
    # lfemm2 = FEMMBase(IntegDomain(loadbdry2, PointRule()))
    # fi2 = ForceIntensity(FFlt[0, 0, -S, 0, 0, 0]);
    bfes = meshboundary(fes)
    el1 = selectelem(fens, bfes, box = [L L d d -Inf Inf], inflate = tolerance)
    loadbdry1 = subset(bfes, el1)
    lfemm1 = FEMMBase(IntegDomain(loadbdry1, GaussRule(1, 2)))
    fi1 = ForceIntensity(FFlt[0, 0, +S/d, 0, 0, 0]);
    el2 = selectelem(fens, bfes, box = [L L -d -d -Inf Inf], inflate = tolerance)
    loadbdry2 = subset(bfes, el2)
    lfemm2 = FEMMBase(IntegDomain(loadbdry2, GaussRule(1, 2)))
    fi2 = ForceIntensity(FFlt[0, 0, -S/d, 0, 0, 0]);
    F = distribloads(lfemm1, geom0, dchi, fi1, 1) + distribloads(lfemm2, geom0, dchi, fi2, 1);

    # @infiltrate
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu =  minimum(dchi.values[:, 3]), maximum(dchi.values[:, 3])
    # @info "Target: $(round.(targetu, digits=8))"

    # Generate a graphical display of displacements and rotations
    scalars = []
    for nc in 1:6
        push!(scalars, ("dchi$nc", deepcopy(dchi.values[:, nc])))
    end
    vectors = []
    push!(vectors, ("U", deepcopy(dchi.values[:, 1:3])))
    push!(vectors, ("UR", deepcopy(dchi.values[:, 4:6])))
    visualize &&  vtkwrite("z_cant-$input-$nrefs-dchi.vtu", fens, fes; scalars = scalars, vectors = vectors)

    # Generate a graphical display of resultants
    nl = selectnode(fens; nearestto = Float64[2.5 1.0 1.0])
    
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :moment, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("m$nc", fld.values))
    end
    pointAstresses = Float64[]
    for nc in 1:3
        fld = fieldfromintegpoints(femm, geom0, dchi, :membrane, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(pointAstresses, fld.values[nl[1], 1]/thickness)
        push!(scalars, ("n$nc", fld.values))
    end
    for nc in 1:2
        fld = fieldfromintegpoints(femm, geom0, dchi, :shear, nc, outputcsys = ocsys)
            # fld = elemfieldfromintegpoints(femm, geom0, dchi, :moment, nc)
        push!(scalars, ("q$nc", fld.values))
    end

    # Visualization
    if visualize
        vtkwrite("z_cant-$input-$nrefs-resultants.vtu", fens, fes; scalars = scalars, vectors = [("u", dchi.values[:, 1:3])])
    end
    return pointAstresses
end

function test_convergence()
    println("--------------------------------------------------------------------------------")
    @info "Z-section cantilever under torsional loading."
    pointAstressX = Float64[]
    for n in [0, 1, 2, 3, 4, 5, ]
        s = _execute("nle5xf3c.inp", n, false)
        push!(pointAstressX, s[1])
    end
    return pointAstressX
end

end # module

using FinEtools
using .LE5_Z_cantilever_examples
@show pointAstressX = LE5_Z_cantilever_examples.test_convergence()
q1, q2, q3 = pointAstressX[end-2:end] ./ phun("MPa")
@show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)