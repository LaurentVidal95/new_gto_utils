using Statistics
using PyCall
using LinearAlgebra
using ProgressMeter
using ThreadsX
using FastGaussQuadrature

################################ GAUSSIAN QUADRATURE INTEGRATION AND SCALAR PRODUCT
function get_gauss_points_and_weights(M)
    # [-1,1]
    ζ, ω = gausslegendre(M)
    # [0,1]
    ζ = (ζ .+ 1) ./2; ω ./= 2;
    
    # Tensorization
    ζ_tensor = []; ω_tensor = [];
    for (i,x) in enumerate(ζ); for (j,y) in enumerate(ζ); for (k,z) in enumerate(ζ);
                push!(ζ_tensor, Float64.([x,y,z]))
                push!(ω_tensor, Float64(ω[i]*ω[j]*ω[k]))
    end; end;  end;
    ω_tensor, ζ_tensor
end

"""
    Convert barycentric coordinates in a given cell to cartesian coordinates
"""
function barycentric_to_cartesian_coords(cell, barycentric_coords)
    # Retrieve cell axes and origin
    x_axis = cell[2] .- cell[1]; y_axis = cell[4] .- cell[1]; z_axis = cell[5] .- cell[1]
    origin = cell[1]
    # Computes the cartesian coordinates
    s,t,u = barycentric_coords
    origin .+ x_axis .*s .+ y_axis .*t .+ z_axis .* u
end

"""
   From the values of a given function on one cell vertices, computes
   the value at point ζ, given in barycentric coordinates of the cell.
"""
function Q1_approximation_in_cell(values_on_vertices, gauss_point)
    e(x,y,z) = [(1-x)*(1-y)*(1-z), x*(1-y)*(1-z), x*y*(1-z), (1-x)*y*(1-z),
                (1-x)*(1-y)*z, x*(1-y)*z, x*y*z, (1-x)*y*z]
    s,t,u = gauss_point
    dot(values_on_vertices, e(s,t,u))
end

function vectors_cubic_transform(cell)
    a0 = cell[1]
    a1 = cell[2] - a0
    a2 = cell[4] - a0
    a3 = cell[5] - a0
    a4 = cell[3] - (a0 + a1 + a2)
    a5 = cell[6] - (a0 + a1 + a3)
    a6 = cell[8] - (a0 + a2 + a3)
    a7 = cell[7] - (a0 + a1 + a2 + a3 + a4 + a5 + a6)

    # set a0 to 0, but does not import for computing jacobian
    0*a0, a1, a2, a3, a4, a5, a6, a7
end

function jacobian_cubic_transform(vectors, ζ)
    a0, a1, a2, a3, a4, a5, a6, a7 = vectors
    x, y, z = ζ
    # J is the Jacobian of the transform from the unit cube to the actual fucked
    # up cell
    J = hcat(a1 + y*a4 + z*a5 + y*z*a7,
             a2 + x*a4 + z*a6 + x*z*a7,
             a3 + x*a5 + y*a6 + x*y*a7)
    abs(det(J)) # |J(x)| for change of variables
end

function volume_cubic_cell(cell)
    x = cell[2] .- cell[1]
    y = cell[4] .- cell[1]
    z = cell[5] .- cell[1]
    det(hcat(x,y,z))
end

"""
   Computes the 3D integral of the Q1 approximation of a function (e.g. Ψ_ex) on one cubic cell, given
   its values on the vertices of the cell
"""
function int3D_on_cubic_cell(cell, Ψ_ex_on_cell, weights, gauss_points)
    # Volume
    volume = volume_cubic_cell(cell)

    int = zero(typeof(Ψ_ex_on_cell[1]))
    for (ω, ζ) in zip(weights, gauss_points)
        Ψ_ζ_Q1 = Q1_approximation_in_cell(Ψ_ex_on_cell, ζ)
        int += ω*Ψ_ζ_Q1
    end
    int*volume
end

"""
   Computes the scalar product of the Q1 approximation of a given function (e.g. Ψ_ex) against
   a function given as an analytic formula, on a given cubic cell.
"""
function scalar_product_on_cubic_cell(cell, Ψ_ex_on_cell, Χμ::Function, weights, gauss_points)
    # Volume
    volume = volume_cubic_cell(cell)

    integral = zero(typeof(Ψ_ex_on_cell[1]))
    for (ω, (s,t,u)) in zip(weights, gauss_points)
        ao_values = Χμ(barycentric_to_cartesian_coords(cell, (s,t,u)))
        Q1_value_of_Ψ_ex = Q1_approximation_in_cell(Ψ_ex_on_cell, (s,t,u))
        integral += ω*Q1_value_of_Ψ_ex*ao_values
    end
    integral*volume
end

function scalar_product_on_cubic_cell(cell, Χν::Function, Χμ::Function, weights, gauss_points)
    # Volume
    volume = volume_cubic_cell(cell)

    integral = zero(Float64)
    for (ω, (s,t,u)) in zip(weights, gauss_points)
        ao_values_1 = Χμ(barycentric_to_cartesian_coords(cell, (s,t,u)))
        ao_values_2 = Χν(barycentric_to_cartesian_coords(cell, (s,t,u)))
        integral += ω*ao_values_1*ao_values_2
    end
    integral*volume
end

function scalar_product_on_cubic_cell(cell, Ψ1, Ψ2::AbstractArray, weights, gauss_points)
    # Volume
    volume = volume_cubic_cell(cell)

    integral = zero(typeof(Ψ1[1][1]))
    for (ω, (s,t,u)) in zip(weights, gauss_points)
        Q1_value_of_Ψ1 = Q1_approximation_in_cell(Ψ1, (s,t,u))
        Q1_value_of_Ψ2 = Q1_approximation_in_cell(Ψ2, (s,t,u))
        integral += ω*Q1_value_of_Ψ1*Q1_value_of_Ψ2
    end
    integral*volume
end


"""
   Same than the preceding functions but on the whole mesh
"""
function int3D_on_cubic_mesh(mesh, Ψ_ex_on_mesh; M=3)
    println(M)
    weights, gauss_points = get_gauss_points_and_weights(M)
    int3D = zero(Float64)
    for (cell, Ψ_ex_on_cell) in zip(mesh, Ψ_ex_on_mesh)
        int3D += int3D_on_cubic_cell(cell, Ψ_ex_on_cell, weights, gauss_points)
    end
    int3D
end

function scalar_product_on_cubic_mesh(mesh, Ψ_ex_on_mesh, Χμ::Function; M=3)
    weights, gauss_points = get_gauss_points_and_weights(M)
    int3D = sum( ThreadsX.map(zip(mesh,Ψ_ex_on_mesh)) do (cell, Ψ_ex_on_cell)
                 scalar_product_on_cubic_cell(cell, Ψ_ex_on_cell, Χμ, weights, gauss_points)
                 end)
    int3D
end

function scalar_product_on_cubic_mesh(mesh, Χν::Function, Χμ::Function; M=3)
    weights, gauss_points = get_gauss_points_and_weights(M)
    int3D = sum( ThreadsX.map(mesh) do (cell)
                 scalar_product_on_cubic_cell(cell, Χν, Χμ, weights, gauss_points)
                 end)
    int3D
end




function scalar_product_on_cubic_mesh(mesh, Ψ1, Ψ2::AbstractArray; M=3)
    weights, gauss_points = get_gauss_points_and_weights(M)
    int3D = sum( ThreadsX.map(zip(mesh,Ψ1,Ψ2)) do (cell, Ψ1_on_cell, Ψ2_on_cell)
                 scalar_product_on_cubic_cell(cell, Ψ1_on_cell, Ψ2_on_cell, weights, gauss_points)
                 end)
    int3D
end

function L2_norm_on_mesh(mesh, f_on_mesh; M=3)
    square_f_on_mesh = ThreadsX.map(f_on_cell-> map(x->norm(x)^2, f_on_cell), f_on_mesh)
    √(int3D_on_cubic_mesh(mesh, square_f_on_mesh,M=M))
end
