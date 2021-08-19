# Julia parser to compute MO in the spatial basis from the MO in given gaussian
# basis from QP2. It needs several input .dat files
# - ao_param.dat :
#   * ao_num = numbers of orbitals
#   * ao_prim_num_max = maximum number of primitive gaussians for AO
#   * nucl_num = number of nuclei
# - nuclei.dat :
#   * nucl_label = label of atoms
#   * nucl_coordinates = coordinates of each nucleus (in Bohrs)
# - mo_coef.dat :
#   contains MO coefficients in a Nb x N matrix, where Nb = ao_num and N is the
#   number of MO computed
# - ao_coef_normalized.dat :
#   contains normalized cj*Nj coefficients in a Nb x m matrix, where Nb = ao_num
#   and m = ao_prim_num_max
# - ao_power.dat :
#   contains power of (x-Rx), (y-Ry), (z-Rz) for each orbitals in a Nb x 3
#   matrix, where Nb = ao_num and 3 is for each coordinate (x,y,z)
# - ao_expo.dat :
#   contains exponants μj in a Nb x m matrix, where Nb = ao_num and
#   m = ao_prim_num_max

using DelimitedFiles

function parse_basis(dir)
    
    println("----------------------------------------")
    println("Loading parameters...\n")

    # Load parameters
    params = readdlm("$(dir)/ao_param.dat", '\n')

    ao_basis = params[1]    # AO basis
    println("Basis: $(ao_basis)")

    Nb = params[2]          # number of AO
    println("Size of basis: $(Nb)")

    m = params[3]    # maximum number of gaussian primitives
    nucl_num = params[4]    # number of nuclei

    # Nuclei label and positions in Bohrs
    nuclei = readdlm("$(dir)/nuclei.dat", ',', '\n')
    nuclei_labels = nuclei[:,1]
    R = nuclei[:,2:4]
    println("Atoms: $(nuclei_labels)")
    println("Positions:")
    display(R)

    # MO coefs
    mo_coefs = readdlm("$(dir)/mo_coef.dat", ',', Float64, '\n')
    println("Size of MO coef matrix:   $(size(mo_coefs))")

    # AO overlap
    ao_overlap = readdlm("$(dir)/ao_overlap.dat", ',', Float64, '\n')

    # AO nuclei index
    ao_nuclei_index = readdlm("$(dir)/ao_nucl_index.dat", ',', Int64, '\n')

    # normalized AO coefs
    ao_coefs = readdlm("$(dir)/ao_coef_normalized.dat", ',', Float64, '\n')
    println("Size of AO coef matrix:   $(size(ao_coefs))")

    # A0 exponants
    σ = readdlm("$(dir)/ao_expo.dat", ',', Float64, '\n')
    println("Size of exponants matrix: $(size(σ))")

    # powers
    p = readdlm("$(dir)/ao_power.dat", ',', Int64, '\n')
    println("Size of powers matrix:    $(size(p))")

    println("\nParameters loaded !")
    println("----------------------------------------")

    # function to parse the μ-th AO from Gaussian primitive basis to space grid at
    # pos = [x,y,z]
    function χ(μ, pos)
        χμ = 0.0
        atom_id = ao_nuclei_index[μ]
        for j = 1:m
            val = ao_coefs[μ, j]
            for (k, r) in enumerate(pos)
                val *= (r - R[atom_id, k])^p[μ, k] * exp(- σ[μ, j] * (r - R[atom_id, k])^2)
            end
            χμ += val
        end
        χμ
    end

    # function to parse the i-th MO from Gaussian basis to space grid at pos=[x,y,z]
    function ϕ(i, pos)
        ϕi = 0.0
        for μ = 1:Nb
            ϕi += mo_coefs[μ, i] * χ(μ, pos)
        end
        ϕi
    end

    χ, ϕ, Nb, ao_overlap
end
