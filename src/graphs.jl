
# a graphplot takes in either an (N x N) adjacency matrix
#   note: you may want to pass node weights to markersize or marker_z
# A graph has N nodes where adj_mat[i,j] is the strength of edge i --> j.  (adj_mat[i,j]==0 implies no edge)


@userplot GraphPlot

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.
function spectral_graph(adjmat::AbstractMatrix)
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
    vec(v[2,:]), vec(v[3,:]), vec(v[4,:])
end


@recipe function f(g::GraphPlot; dim = 2, T = Float64)
    @assert dim in (2, 3)
    delete!(d, :dim)
    delete!(d, :T)

    adjmat = g.args[1]
    n, m = size(adjmat)
    x, y, z = spectral_graph(adjmat)

    # create a series for the line segments
    if get(d, :linewidth, 1) > 0
        @series begin
            # skipped when user overrides linewidth to 0
            # we want to build new lx/ly/lz for the lines
            # note: we only do the lower triangle

            lx, ly, lz = zeros(T,0), zeros(T,0), zeros(T,0)
            line_z = zeros(T,0)
            for i=1:n, j=i+1:n
                aij = adjmat[i,j]
                if aij ≉ 0
                    append!(lx, T[x[i], x[j], NaN])
                    append!(ly, T[y[i], y[j], NaN])
                    if dim == 3
                        append!(lz, T[z[i], z[j], NaN])
                    end
                    append!(line_z, T[aij, aij])
                    # TODO: when supported, add line width for this line segment
                end
            end
            series_annotations := []
            # linecolor --> :black
            linewidth --> 1
            line_z := line_z
            markershape := :none
            markercolor := :black
            primary := false
            dim==3 ? (lx, ly, lz) : (lx, ly)
        end
    end

    seriestype := (dim==3 ? :scatter3d : :scatter)
    linewidth := 0
    linealpha := 0
    foreground_color_border := nothing
    grid := false
    legend := false
    ticks := nothing
    # markersize := node_weight
    dim==3 ? (x, y, z) : (x, y)
end
