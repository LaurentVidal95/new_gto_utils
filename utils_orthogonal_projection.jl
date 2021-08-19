using LinearAlgebra


###################### With d(C) = 1 - C'*S*C
"""
    Computes the distance with the formula d = 1 - V'inv(S)V
"""
function scalar_products_with_AO_basis(mesh, Ψ_ex, Χμs, Nb; M=3)
    [scalar_product_on_cubic_mesh(mesh, Ψ_ex, X->Χμs(i,X), M=M) for
     i in 1:Nb]
end

distance_to_AO_basis(C,S) = 1 - C'*S*C

function distance_to_AO_basis(mesh, Ψ_ex, Χμs, Nb, S; M=3)
    V = scalar_products_with_AO_basis(mesh, Ψ_ex, Χμs, Nb, M=M)
    C = (S\V)
    distance_to_AO_basis(C,S)
end


###################### With approx ||Ψ_ex - ∑cμΧμ||
"""
    Compute the distance with ||Ψ_ex - ∑cμΧμ||_{L^2}
"""
Ψ_proj(X, C, Χμs, Nb) = sum([C[μ]*Χμs(μ,X) for μ in 1:Nb])

function distance_to_AO_basis_with_approx_L2_norm(Ψ_ex, Ψ_proj_on_mesh, mesh;M=3)
    # Correct eventual phase opposition
    Ψ_ex = test_phase_opposition(Ψ_ex, Ψ_proj_on_mesh) .* Ψ_ex
    Ψ_ex_minus_proj_on_mesh = [values_ex .- values_projs for (values_ex, values_projs) in
                               zip(Ψ_ex,Ψ_proj_on_mesh)]
    L2_norm_on_mesh(mesh, Ψ_ex_minus_proj_on_mesh, M=M)
end

function test_phase_opposition(Ψ1, Ψ2)
    test = sign(Ψ1[end][1]*Ψ2[end][1])
    # If the tested values are zero, code a better function
    (test == 0) && (@error("Error testing phase shift: null values"))
    test
end
