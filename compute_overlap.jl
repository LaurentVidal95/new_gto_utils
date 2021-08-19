# Dependancies
include("include_dependencies.jl")

# Extract Ψ_ex data
DFTFE_pvtu_file = "./H2_ion/waveFunctionOutputFolder/wfcOutput_master.pvtu"
mesh, N_cells, Ψi_dftfe = extract_mesh_data(DFTFE_pvtu_file)

Ψ0_dftfe = Ψi_dftfe[1] # Ground state

# List all basis
basis_dir = "basis_data/" # "basis_data/"
basis_list = readdir(basis_dir)

overlap_matrices_qp = []
overlap_matrices = []

# Compute the distancess for all basis
for basis in basis_list
    # Extract analytc
    Χμ_list, φk_list, Nb, S = parse_basis("$(basis_dir)/$(basis)")
    push!(overlap_matrices_qp, S)

    # compute overlap matrix
    S2 = zeros(Nb, Nb)
    for i = 1:Nb
        for j = 1:Nb
            println((i,j))
            S2[i,j] = scalar_product_on_cubic_mesh(mesh,
                                                   X->Χμ_list(i,X),
                                                   X->Χμ_list(j,X), M=3)
        end
    end
    push!(overlap_matrices, S2)
    STOP
end

