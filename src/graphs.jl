
# I'm giving this a silly name because eventually I want to use Graphs.jl and I don't want clashes
type MyGraph{N <: AbstractVector, E <: AbstractMatrix}
    nodewgt::N
    adjmat::E
end

# create an uninitialized cell then add the args
function Base.cell(arg1::AbstractVector, args...)
    c = cell(length(args) + 1)
    c[1] = arg1
    for i in 1:length(args)
        c[i+1] = args[i]
    end
    c
end

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf
@recipe function f(g::MyGraph, dim::Integer = 2)
    @assert dim in (2, 3)

    A = g.adjmat
    n, m = size(A)
    @assert n == m

    # D is a diagonal matrix with the degrees (total weight for that node) on the diagonal
    deg = vec(sum(A,1))
    D = diagm(deg)

    # Laplacian (L = D - A)
    L = Float64[i == j ? deg[i] : -A[i,j] for i=1:n,j=1:n]

    # get the matrix of eigenvectors
    v = eig(L, D)[2]

    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    x, y, z = vec(v[2,:]), vec(v[3,:]), vec(v[4,:])

    :linewidth --> 1
    if d[:linewidth] > 0
        # skipped when user overrides linewidth to 0
        # we want to build new lx/ly/lz for the lines
        lx, ly, lz = zeros(0), zeros(0), zeros(0)
        for i=1:n, j=1:n
            if A[i,j] ≉ 0
                append!(lx, Float64[x[i], x[j], NaN])
                append!(ly, Float64[y[i], y[j], NaN])
                append!(lz, Float64[z[i], z[j], NaN])
                # TODO: when supported, add line color/width for this line segment
            end
        end
    end

    :markershape --> :circle
    :markersize --> 10 * g.nodewgt

    if dim == 2
        :linetype --> [:scatter :path]
        cell(x, lx), cell(y, ly)
    else
        :linetype --> [:scatter3d :path]
        cell(x lx), cell(y ly), cell(z, lz)
    end
end
