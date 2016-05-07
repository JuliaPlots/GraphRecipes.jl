
# I'm giving this a silly name because eventually I want to use Graphs.jl and I don't want clashes
type MyGraph{N <: AbstractVector, E <: AbstractMatrix}
    nodewgt::N
    adjmat::E
end

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
@recipe function f(g::MyGraph)

    A = g.adjmat
    n, m = size(A)
    @assert n == m
    deg = vec(sum(A,1))

    # normalized Laplacian (N = I - D^(-1/2) * A * D^(-1/2))
    # TODO check for division by zero
    normlap = Float64[(i == j ? 1 : 0) - (A[i,j] / sqrt(deg[i]*deg[j])) for i=1:n, j=1:n]
    @show normlap

    rand(n)
end
