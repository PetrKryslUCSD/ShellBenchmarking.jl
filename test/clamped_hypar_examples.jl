"""
Clamped hyperbolic paraboloid under self weight.

From: CMES, vol.49, no.2, pp.81-110, 2009

The problem considered in this section is that of a hyperbolic paraboloid shell,
clamped along one side and free on three edges and loaded by self-weight
(Figure 16). This is a pure bending dominated problem and known to be a very
hard test for locking behaviour as suggested in References (Chapelle and Bathe,
1998; Bathe, Iosilevich, and Chapelle, 2000). 

The shell geometry is described by the
equation: z = x^2 −y^2 ; (x,y) ∈ − L/2 ; L/2

Results:

Table 2 in Bathe, Iosilevich, and Chapelle, 2000 (computed 
with one fixed mesh using MITC16 elements: these are not asymptotic results).
t/L   Strain energy  Displacement
1/100 1.6790e-3 9.3355e-5
1/1000 1.1013e-2 6.3941e-3
1/10,000 8.9867e-2 5.2988e-1

See also the table in  
An Improved Quadrilateral Flat Element with Drilling
Degrees of Freedom for Shell Structural Analysis
H. Nguyen-Van1 , N. Mai-Duy1 and T. Tran-Cong1
CMES, vol.49, no.2, pp.81-110, 2009.

"""
module clamped_hypar_examples

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute_half(orientation = :a, tL_ratio = 1/100, n = 32, visualize = false)
    # Parameters:
    E = 2.0e11 * phun("Pa");
    nu = 0.3;
    L = 1.0 * phun("m")
    g = 8000 * phun("N/m^3")
    thickness = tL_ratio * L
    formul = FEMMShellT3FFModule

    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L/2,n,Int(round(n/2)),orientation);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # femm.drilling_stiffness_scale = 0.1
    # femm.mult_el_size = 0.2
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    # Symmetry Plane
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1), thickness))
    fi = ForceIntensity(Float64[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    solve_blocked!(dchi, K, F)
    targetu = dchi.values[nl, 3][1]
    U = gathersysvec(dchi, DOF_KIND_ALL)
    targetse = 1/2*U'*F*2
    @info "Deflection at A: $(round(targetu, digits = 9))"
    @info "Strain Energy (entire structure): $(targetse) "
        # Visualization
    if visualize

        vtkwrite("clamped_hypar-$(orientation)-$(n).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu, targetse
end

function test_convergence(orientation = :a)
    println("--------------------------------------------------------------------------------")
    @info "Clamped hyperbolic paraboloid under self weight."
    tL_ratios = [1/100, 1/1000, 1/10000]; 
    
    ns = [4, 8, 16, 32, 64, 128, 256, 512]
    all_uresults = []
    all_seresults = []
    for (tL_ratio, ) in zip(tL_ratios, )
        @info "Clamped hypar, t/L=$(tL_ratio)"
        uresults = Float64[]
        seresults = Float64[]
        for n in ns
            u, se = _execute_half(orientation, tL_ratio, n, false)
            push!(uresults, u)
            push!(seresults, se)
        end   
        push!(all_uresults, uresults)
        push!(all_seresults, seresults)
    end

    return ns, all_uresults, all_seresults
end

end # module

using .clamped_hypar_examples

let
    for orientation in (:a, :b)
        ns, all_uresults, all_seresults = clamped_hypar_examples.test_convergence(orientation)
    
        @show ns
        @show all_uresults
        @show all_seresults
        for i in 1:3
            q1, q2, q3 = all_uresults[i][end-2:end]
            @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
            q1, q2, q3 = all_seresults[i][end-2:end]
            @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
        end
    end 
end 
