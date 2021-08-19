function test_cube(cell)
    x_axis = cell[2] .- cell[1]; y_axis = cell[4] .- cell[1]; z_axis = cell[5] .- cell[1]
    abs(dot(x_axis,y_axis))/(norm(x_axis)*norm(y_axis)),abs(dot(x_axis,z_axis))/(norm(x_axis)*norm(z_axis)), abs(dot(y_axis,z_axis))/(norm(y_axis)*norm(z_axis))
end

function compute_volume_hexa(cell)
    x1, x2, x3, x4, x5, x6, x7, x8 = cell
    v1 = hcat(x7-x1, (x2+x8)-(x5+x6), -(x2-x8)+(x5-x6))
    v2 = hcat(x7-x1, (x2+x8)-(x3+x4), (x2-x8)+(x3-x4))
    (det(v1) + det(v2))/12
end

function compute_volume_hexa_bis(cell)
    x1, x2, x3, x4, x5, x6, x7, x8 = cell
    v1 = hcat(x7-x1, x2-x1, x3-x6)
    v2 = hcat(x7-x1, x5-x1, x6-x8)
    v3 = hcat(x7-x1, x4-x1, x8-x3)
    (det(v1) + det(v2) + det(v3))/6
end

function test_mesh(mesh; out=false)
    cubic = true
    for (i,cell) in enumerate(mesh)
        if sum(test_cube(cell)) > eps(typeof(mesh[1][1][1]))
            cubic = false
            if out
                println("cellule : $i, valeurs : $(sum(test_cube(cell)))")
                println("$(cell[1][1])")
            end
        end
    end
    println(cubic)
    return cubic
end
