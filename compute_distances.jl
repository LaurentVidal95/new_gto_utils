# Dependancies
include("include_dependencies.jl")

# Extract Ψ_ex data
DFTFE_pvtu_file = "./H2_DFTFE_solution/low_precision_data/wfc_H2_Q1_base_mesh_10_atom_mesh_015/wfcOutput_master.pvtu"
mesh, N_cells, Ψi_dftfe = extract_mesh_data(DFTFE_pvtu_file)

Ψ0_dftfe = Ψi_dftfe[1] # Ground state

# List all basis
basis_dir = "basis_data/" # "basis_data/"
basis_list = readdir(basis_dir)

list_distances = Float64[]

# Compute the distancess for all basis
for basis in basis_list
    # Extract analytc
    Χμ_list, φk_list, Nb, S = parse_basis("$(basis_dir)/$(basis)")

    # Compute Ψ0_dftfe projected on the AO basis
    push!(list_distances, distance_to_AO_basis(mesh, Ψ0_dftfe, Χμ_list, Nb, S))
end

f = open("list_distances.dat","w")
for d in list_distances
    write(f,"$d ")
end
close(f)
