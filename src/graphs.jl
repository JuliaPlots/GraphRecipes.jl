
# # I'm giving this a silly name because eventually I want to use Graphs.jl and I don't want clashes
# type MyGraph{N <: AbstractVector, E <: AbstractMatrix}
#     nodewgt::N
#     adjmat::E
# end

@userplot GraphPlot

# # TODO: add this to RecipesBase?  Should probably just name it something else, even though
# #       it shares the same goal (construction of a Vector{Any})
# # create an uninitialized cell then add the args
# function Base.cell(arg1::AbstractVector, args...)
#     c = cell(length(args) + 1)
#     c[1] = arg1
#     for i in 1:length(args)
#         c[i+1] = args[i]
#     end
#     c
# end

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.
@recipe function f(g::GraphPlot; dim = 2)
    @assert dim in (2, 3)
    delete!(d, :dim)
    # is3d = dim == 3

    # the args are expected to be either adjmat or (node_weight, adjmat)
    node_weight, adjmat = if length(g.args) == 1
        adjmat = g.args[1]
        zeros(size(adjmat,1)), adjmat
    else
        g.args
    end

    n, m = size(adjmat)
    @assert n == m

    # D is a diagonal matrix with the degrees (total weight for that node) on the diagonal
    deg = vec(sum(adjmat,1))
    D = diagm(deg)

    # Laplacian (L = D - adjmat)
    L = Float64[i == j ? deg[i] : -adjmat[i,j] for i=1:n,j=1:n]

    # get the matrix of eigenvectors
    v = eig(L, D)[2]

    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    x, y, z = vec(v[2,:]), vec(v[3,:]), vec(v[4,:])

    # create a series for the line segments
    if get(d, :linewidth, 1) > 0
        @series begin
            # skipped when user overrides linewidth to 0
            # we want to build new lx/ly/lz for the lines
            lx, ly, lz = zeros(0), zeros(0), zeros(0)
            for i=1:n, j=1:n
                if adjmat[i,j] ≉ 0
                    append!(lx, Float64[x[i], x[j], NaN])
                    append!(ly, Float64[y[i], y[j], NaN])
                    append!(lz, Float64[z[i], z[j], NaN])
                    # TODO: when supported, add line color/width for this line segment
                end
            end
            series_annotations := []
            linecolor --> :black
            linewidth --> 1
            markershape := :none
            markercolor := :black
            label := ""
            primary := false
            dim==3 ? (lx, ly, lz) : (lx, ly)
        end
    end

    seriestype := (dim==3 ? :scatter3d : :scatter)
    linewidth := 0
    linealpha := 0
    dim==3 ? (x, y, z) : (x, y)

    # linewidth --> 1
    # markershape --> [:circle :none]

    # # TODO: adjust marker size to nodewgt
    # # :markersize --> cell(10 * g.nodewgt, 0)

    # if d[:linewidth] > 0
    #     # skipped when user overrides linewidth to 0
    #     # we want to build new lx/ly/lz for the lines
    #     lx, ly, lz = zeros(0), zeros(0), zeros(0)
    #     for i=1:n, j=1:n
    #         if A[i,j] ≉ 0
    #             append!(lx, Float64[x[i], x[j], NaN])
    #             append!(ly, Float64[y[i], y[j], NaN])
    #             append!(lz, Float64[z[i], z[j], NaN])
    #             # TODO: when supported, add line color/width for this line segment
    #         end
    #     end

    #     d[:linewidth] = [0 d[:linewidth]]
    #     :linetype --> (is3d ? [:scatter3d :path3d] : [:scatter :path])
    #     is3d ? (cell(x, lx), cell(y, ly), cell(z, lz)) : (cell(x, lx), cell(y, ly))
    # else
    #     # no need for line segments
    #     :linetype --> (is3d ? :scatter3d : :scatter)
    #     is3d ? (x, y, z) : (x, y)
    # end
end
