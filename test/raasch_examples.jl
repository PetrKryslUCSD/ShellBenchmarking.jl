"""
The geometry consists of a “hook” in the form of a curved strip rigidly clamped
at one end and loaded with a unit in-plane shear along the width at the other
end. It has two circular segments that are connected at the tangent point. The
smaller segment has a mean radius of 0.3556 m (14 inches) and spans 60° from
the clamped end to the tangent point. The larger segment spans 150° from the
tangent point to the free end and has a mean radius of 1.1684 m (46 inches).
The hook is 0.0508 m (2 inches) thick and 0.508 m (20 inches) wide, modeled as
linear elastic with an elastic modulus of 22.77 MPa (3300 psi) and a Poisson's
ratio of 0.35. In most tests the shear force is applied through the use of a
distributing coupling constraint. The coupling constraint provides coupling
between a reference node on which the load is prescribed and the nodes located
on the free end. The distributed nodal loads on the free end are equivalent to
a uniformly distributed load of 8.7563 N/m (0.05 lb/in). In two of the tests an
equivalent shear force is applied as a distributed shear traction instead.


"""
module raasch_examples

using LinearAlgebra, Statistics
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellT3FFModule: num_normals
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

bench_sol = 5.022012648671993*phun("in"); # 20-node hex results

function _execute(nL = 9, nW = 1, visualize = true)
    E = 3300.0 * phun("psi");
    nu = 0.35;
    thickness  =  2.0 * phun("in");
    tolerance = thickness/20
    R1 = 14.0*phun("in");
    R2 = 46.0*phun("in");
    q = 0.05 * phun("lbf/in")
    formul = FEMMShellT3FFModule
    
    fens, fes = T3block(210.0, 20.0 * phun("in"), nL, nW)
    fens.xyz = xyz3(fens)
    for k in 1:count(fens)
        a, z = fens.xyz[k, :]
        if a <= 60.0
            fens.xyz[k, :] .= (R1*sin(a*pi/180), R1*(1-cos(a*pi/180)), z)
        else
            a = a - 60.0 + 30.0
            xc = R1*sin(60*pi/180) + R2*cos(30*pi/180)
            yc = R1*(1-cos(60*pi/180)) - R2*sin(30*pi/180)
            fens.xyz[k, :] .= (xc-R2*cos(a*pi/180), yc+R2*sin(a*pi/180), z)
        end
    end
        
    @show count(fens), count(fes)

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
    # Clamped end
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, box = [97.96152422706632 97.96152422706632 -16 -16 0 20].*phun("in"), inflate = 1.0e-6)
    lfemm = FEMMBase(IntegDomain(subset(bfes, l1), GaussRule(1, 2)))
    fi = ForceIntensity(Float64[0, 0, q, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    @assert  isapprox(sum(F)/phun("lbf"), 1.0)
    
    # Solve
    solve_blocked!(dchi, K, F)
    nl = selectnode(fens; box = Float64[97.96152422706632 97.96152422706632 -16 -16 0 20].*phun("in"), inflate = 1.0e-6)
    targetu =  mean(dchi.values[nl, 3])
    @info "Solution: $(round(targetu, digits=8)/phun("in")) [in],  $(round(targetu/bench_sol, digits = 4)*100)%"

    # Visualization
    if visualize
        vtkwrite("raasch-geometry.vtu", fens, fes)
        scattersysvec!(dchi, (R1/2)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -R1]; [R1 R1 R1]]),
            plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu
end

function test_convergence()
    println("--------------------------------------------------------------------------------")
    @info "Raasch hook."
    results = []
    for n in 1:7
        push!(results, _execute(9*2^(n-1), 1*2^(n-1), false))
    end
    return results
end

end # module

using FinEtools
using .raasch_examples
results = raasch_examples.test_convergence()
q1, q2, q3 = results[end-2:end] ./ phun("in")
@show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
