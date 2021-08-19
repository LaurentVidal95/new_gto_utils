function extract_mesh_data(file; N=1, raw = false)
    # Check presence of the file and read data as python object
    @assert(isfile(file))
    vs = pyimport("pyvista")
    mesh = vs.read(file)
    
    # Extract raw data as julia arrays
    cells, points = mesh.cells, mesh.points
    N_cells = Int64(length(cells)/9)

    (raw) && (return cells, N_cells, points)
    
    # Reshape so that cells = [[cell_1_coords], [cell_2_coords], etc.. ]
    cells_indexes = [cells[9*(i-1)+2:9*(i-1)+9] .+ 1 for i in 1:N_cells]
    cells_coords = [[points[id,:] for id in indexes] for indexes in cells_indexes]

    # Extract values as [Ψ1, Ψ2, ...] where Ψi = [[values of Ψi at cell 1], etc.. ]
    Ψs = mesh.point_arrays; N_Ψs = length(Ψs);
    values = []
    for i in 1:min(N_Ψs, N)
        Ψi = Ψs[i]
        tmp = [[Ψi[id] for id in indexes] for indexes in cells_indexes]
        push!(values, tmp)
    end
    
    cells_coords, N_cells, values
end
