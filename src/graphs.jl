
# -----------------------------------------------------

function get_source_destiny_weight{T}(mat::AbstractArray{T,2})
    nrow, ncol = size(mat)    # rows are sources and columns are destinies

    nosymmetric = !issymmetric(mat) # plots only triu for symmetric matrices
    nosparse = !issparse(mat) # doesn't plot zeros from a sparse matrix

    L = length(mat)

    source  = Array(Int, L)
    destiny = Array(Int, L)
    weights  = Array(T, L)

    idx = 1
    for i in 1:nrow, j in 1:ncol
        value = mat[i, j]
        if !isnan(value) && ( nosparse || value != zero(T) ) # TODO: deal with Nullable

            if i < j
                source[idx]  = i
                destiny[idx] = j
                weights[idx]  = value
                idx += 1
            elseif nosymmetric && (i > j)
                source[idx]  = i
                destiny[idx] = j
                weights[idx]  = value
                idx += 1
            end

        end
    end

    resize!(source, idx-1), resize!(destiny, idx-1), resize!(weights, idx-1)
end

function get_source_destiny_weight(source::AVec, destiny::AVec)
    if length(source) != length(destiny)
        throw(ArgumentError("Source and destiny must have the same length."))
    end
    source, destiny, Float64[ 1.0 for i in source ]
end

function get_source_destiny_weight(source::AVec, destiny::AVec, weights::AVec)
    if !(length(source) == length(destiny) == length(weights))
        throw(ArgumentError("Source, destiny and weights must have the same length."))
    end
    source, destiny, weights
end


# -----------------------------------------------------


function get_adjacency_matrix(source::AVec{Int}, destiny::AVec{Int}, weights::AVec)
    n = max(maximum(source), maximum(destiny))
    full(sparse(source, destiny, weights, n, n))
end

function make_symmetric(A::AMat)
    A = copy(A)
    for i=1:size(A,1), j=i+1:size(A,2)
        A[i,j] = A[j,i] = A[i,j]+A[j,i]
    end
    A
end

function compute_laplacian(adjmat::AMat, node_weights::AVec)
    n, m = size(adjmat)
    # @show size(adjmat), size(node_weights)
    @assert n == m == length(node_weights)

    # scale the edge values by the product of node_weights, so that "heavier" nodes also form
    # stronger connections
    adjmat = adjmat .* sqrt(node_weights * node_weights')

    # D is a diagonal matrix with the degrees (total weights for that node) on the diagonal
    deg = vec(sum(adjmat,1)) - diag(adjmat)
    D = diagm(deg)

    # Laplacian (L = D - adjmat)
    L = Float64[i == j ? deg[i] : -adjmat[i,j] for i=1:n,j=1:n]

    L, D
end


# -----------------------------------------------------
# -----------------------------------------------------

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.
function spectral_graph(adjmat::AMat; node_weights::AVec = ones(size(adjmat,1)), kw...)
    adjmat = make_symmetric(adjmat)
    L, D = compute_laplacian(adjmat, node_weights)

    # get the matrix of eigenvectors
    v = eig(L, D)[2]

    # x, y, and z are the 2nd through 4th eigenvectors of the solution to the
    # generalized eigenvalue problem Lv = λDv
    vec(v[2,:]), vec(v[3,:]), vec(v[4,:])
end

function spectral_graph(source::AVec{Int}, destiny::AVec{Int}, weights::AVec; kw...)
    spectral_graph(get_adjacency_matrix(source, destiny, weights); kw...)
end


# -----------------------------------------------------

# Axis-by-Axis Stress Minimization -- Yehuda Koren and David Harel
# See: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.437.3177&rep=rep1&type=pdf


# # NOTES:
# #   - dᵢⱼ = the "graph-theoretical distance between nodes i and j"
# #         = Aᵢⱼ
# #   - kᵢⱼ = dᵢⱼ⁻²
# #   - b̃ᵢ = ∑ᵢ≠ⱼ ((x̃ⱼ ≤ x̃ᵢ ? 1 : -1) / dᵢⱼ)
# #   - need to solve for x each iteration: Lx = b̃

# # Solve for one axis at a time while holding the others constant.
# # dims is 2 (2D) or 3 (3D).  free_dims is a vector of the dimensions to update (for example if you fix y and solve for x)
# function by_axis_stress_graph(adjmat::AMat, node_weights::AVec = ones(size(adjmat,1));
#                               dims = 2, free_dims = 1:dims,
#                               x = rand(length(node_weights)),
#                               y = rand(length(node_weights)),
#                               z = rand(length(node_weights)))
#     adjmat = make_symmetric(adjmat)
#     L, D = compute_laplacian(adjmat, node_weights)

#     n = length(node_weights)
#     maxiter = 100 # TODO: something else

#     @assert dims == 2

#     @show adjmat L

#     for _ in 1:maxiter
#         x̃ = x
#         b̃ = Float64[sum(Float64[(i==j || adjmat[i,j] == 0) ? 0.0 : ((x̃[j] <= x̃[i] ? 1.0 : -1.0) / adjmat[i,j]) for j=1:n]) for i=1:n]
#         @show x̃ b̃
#         x = L \ b̃

#         xdiff = x - x̃
#         @show norm(xdiff)
#         if norm(xdiff) < 1e-4
#             info("converged. norm(xdiff) = $(norm(xdiff))")
#             break
#         end
#     end
#     @show x y
#     x, y, z
# end

norm_ij(X, i, j) = sqrt(sum(Float64[(v[i]-v[j])^2 for v in X]))
stress(X, dist, w, i, j) = w[i,j] * (norm_ij(X, i, j) - dist[i,j])^2
function stress(X, dist, w)
    tot = 0.0
    for i=1:size(X,1), j=1:i-1
        tot += stress(X, dist, w, i, j)
    end
    tot
end


# follows section 2.3 from http://link.springer.com/chapter/10.1007%2F978-3-540-31843-9_25#page-1
# Localized optimization, updates: x
function by_axis_local_stress_graph(adjmat::AMat;
                              node_weights::AVec = ones(size(adjmat,1)),
                              dim = 2, free_dims = 1:dim,
                              x = rand(length(node_weights)),
                              y = rand(length(node_weights)),
                              z = rand(length(node_weights)),
                              maxiter = 300,
                              kw...)
    adjmat = make_symmetric(adjmat)
    n = length(node_weights)

    # graph-theoretical distance between node i and j (i.e. shortest path distance)
    # TODO: calculate a real distance
    dist = estimate_distance(adjmat)
    # @show dist

    # also known as kᵢⱼ in "axis-by-axis stress minimization".  the -2 could also be 0 or -1?
    w = dist .^ -2

    # in each iteration, we update one dimension/node at a time, reducing the total stress with each update
    X = dim == 2 ? (x, y) : (x, y, z)
    laststress = stress(X, dist, w)
    for k in 1:maxiter
        for p in free_dims
            for i=1:n
                numer, denom = 0.0, 0.0
                for j=1:n
                    i==j && continue
                    numer += w[i,j] * (X[p][j] + dist[i,j] * (X[p][i] - X[p][j]) / norm_ij(X, i, j))
                    denom += w[i,j]
                end
                if denom != 0
                    X[p][i] = numer / denom
                end
            end
        end

        # check for convergence of the total stress
        thisstress = stress(X, dist, w)
        if abs(thisstress - laststress) / abs(laststress) < 1e-5
            # info("converged. numiter=$k last=$laststress this=$thisstress")
            break
        end
        laststress = thisstress
    end

    dim == 2 ? (X..., nothing) : X
end

function by_axis_local_stress_graph(source::AVec{Int}, destiny::AVec{Int}, weights::AVec; kw...)
    by_axis_local_stress_graph(get_adjacency_matrix(source, destiny, weights); kw...)
end

# -----------------------------------------------------

if Plots.is_installed("LightGraphs")
    @eval begin
        import LightGraphs

        # TODO: so much wasteful conversion... do better
        function estimate_distance(adjmat::AMat)
            source, destiny, weights = get_source_destiny_weight(sparse(adjmat))
            n = size(adjmat,1)
            g = LightGraphs.Graph(n)
            for (si,di,wi) in zip(source,destiny,weights)
                LightGraphs.add_edge!(g, si, di)
            end
            convert(Matrix{Float64}, hcat(map(i->LightGraphs.dijkstra_shortest_paths(g, i).dists, 1:n)...))
        end
    end
else
    @eval function estimate_distance(adjmat::AMat)
        warn("Install LightGraphs for the best layout calculations.")
        map(a -> a==0 ? 5.0 : a, adjmat)
    end
end

function tree_graph(adjmat::AMat; kw...)
    tree_graph(get_source_destiny_weight(adjmat)...; kw...)
end

function tree_graph(source::AVec{Int}, destiny::AVec{Int}, weights::AVec;
                    node_weights::AVec = ones(max(maximum(source), maximum(destiny))),
                    root::Symbol = :top,  # flow of tree: left, right, top, bottom
                    layers_scalar = 1.0,
                    layers = nothing,
                    dim = 2,
                    kw...)
    extrakw = KW(kw)
    n = length(node_weights)

    # @show root

    # TODO: compute layers, which get bigger as you go away from the root
    if layers == nothing
        # layers = rand(1:4, n)
        layers = compute_tree_layers2(source, destiny, n)
    end

    # reverse direction?
    if root in (:top, :right)
        layers = -layers
    end

    # add noise
    layers = layers + 0.6rand(size(layers)...)

    # TODO: normalize layers somehow so it's in line with distances
    layers .*= layers_scalar
    if dim == 2
        if root in (:top, :bottom)
            extrakw[:y] = layers
            extrakw[:free_dims] = [1]
        elseif root in (:left, :right)
            extrakw[:x] = layers
            extrakw[:free_dims] = [2]
        else
            error("unknown root: $root")
        end
    else
        error("3d not supported")
    end

    # now that we've fixed one dimension, let the stress algo solve for the other(s)
    by_axis_local_stress_graph(get_adjacency_matrix(source, destiny, weights);
                               node_weights=node_weights,
                               dim=dim,
                               extrakw...)
end


function adjlist_and_degrees(source, destiny, n)
    # build a list of children (adjacency list)
    alist = Vector{Int}[Int[] for i=1:n]
    indeg, outdeg = zeros(Int, n), zeros(Int, n)
    for (si,di) in zip(source, destiny)
        push!(alist[si], di)
        indeg[di] += 1
        outdeg[si] += 1
    end
    alist, indeg, outdeg
end

function compute_tree_layers(source, destiny, n)
    alist, indeg, outdeg = adjlist_and_degrees(source, destiny, n)

    # choose root to be the node with lots going out, but few coming in
    netdeg = outdeg - 50indeg
    idxs = sortperm(netdeg, rev=true)
    # rootidx = findmax(netdeg)
    # @show outdeg indeg netdeg idxs alist
    placed = Int[]

    layers = zeros(n)
    for i=1:n
        idx = shift!(idxs)

        # first, place this after its parents
        for j in placed
            if idx in alist[j]
                layers[idx] = max(layers[idx], layers[j] + 1)
            end
        end
        
        # next, shift its children lower
        for j in idxs
            if j in alist[idx]
                layers[j] = max(layers[j], layers[idx] + 1)
            end
        end
        
        push!(placed, idx)
    end
    layers
end

# an alternative algo to pick tree layers... generate a list of roots,
# and for each root, make a pass through the tree (without recurrency)
# and push the children below their parents
function compute_tree_layers2(source, destiny, n)
    alist, indeg, outdeg = adjlist_and_degrees(source, destiny, n)
    roots = filter(i->indeg[i]==0, 1:n)
    if isempty(roots)
        roots = [1]
    end

    layers = zeros(Int,n)
    for i in roots
        shift_children!(layers, alist, Int[], i)
    end

    # now that we've shifted children out, move parents closer to their closest children
    while true
        shifted = false
        for parent=1:n
            if !(isempty(alist[parent]))
                minidx = minimum(layers[child] for child in alist[parent])
                if layers[parent] < minidx - 1
                    shifted = true
                    layers[parent] = minidx - 1
                end
            end
        end
        shifted || break
    end

    layers
end

function shift_children!(layers, alist, placed, parent)
    for idx in alist[parent]
        if !(idx in placed) && layers[idx] <= layers[parent]
            layers[idx] = layers[parent] + 1
        end
    end
    for idx in alist[parent]
        if idx != parent && !(idx in placed)
            push!(placed, idx)
            shift_children!(layers, alist, placed, idx)
        end
    end
end


# -----------------------------------------------------

# TODO: maybe also implement Catmull-Rom Splines? http://www.mvps.org/directx/articles/catmull/

# -----------------------------------------------------

# we want to randomly pick a point to be the center control point of a bezier
# curve, which is both equidistant between the endpoints and normally distributed
# around the midpoint
function random_control_point(xi, xj, yi, yj, curvature_scalar)
    xmid = 0.5 * (xi+xj)
    ymid = 0.5 * (yi+yj)

    # get the angle of y relative to x
    theta = atan((yj-yi) / (xj-xi)) + 0.5pi

    # calc random shift relative to dist between x and y
    dist = sqrt((xj-xi)^2 + (yj-yi)^2)
    dist_from_mid = curvature_scalar * (rand()-0.5) * dist

    # now we have polar coords, we can compute the position, adding to the midpoint
    (xmid + dist_from_mid * cos(theta),
     ymid + dist_from_mid * sin(theta))
end

# -----------------------------------------------------

const _graph_funcs = KW(
    :spectral => spectral_graph,
    :stress => by_axis_local_stress_graph,
    :tree => tree_graph,
)

# a graphplot takes in either an (N x N) adjacency matrix
#   note: you may want to pass node weights to markersize or marker_z
# A graph has N nodes where adj_mat[i,j] is the strength of edge i --> j.  (adj_mat[i,j]==0 implies no edge)

# NOTE: this is for undirected graphs... adjmat should be symmetric and non-negative

@userplot GraphPlot

@recipe function f(g::GraphPlot;
                   dim = 2,
                   free_dims = 1:dim,
                   T = Float64,
                   curves = true,
                   curvature_scalar = 0.2,
                   root = :top,
                   node_weights = nothing,
                   names = [],
                   x = nothing,
                   y = nothing,
                   z = nothing,
                   func = spectral_graph,
                   shorten = 0.1
                  )
    @assert dim in (2, 3)
    _3d = dim == 3

    source, destiny, weights = get_source_destiny_weight(g.args...)
    if !(eltype(source) <: Integer)
        names = unique(sort(vcat(source,destiny)))
        source = Int[findfirst(names, si) for si in source]
        destiny = Int[findfirst(names, di) for di in destiny]
    end
    n = max(maximum(source), maximum(destiny))

    if node_weights == nothing
        node_weights = ones(n)
    end
    @assert length(node_weights) == n

    # do we want to compute coordinates?
    if (_3d && (x == nothing || y == nothing || z == nothing)) || (!_3d && (x == nothing || y == nothing))
        if isa(func, Symbol)
            func = _graph_funcs[func]
        end
        x, y, z = func(
            source, destiny, weights;
            node_weights = node_weights,
            dim = dim,
            free_dims = free_dims,
            root = root
        )
    end

    # create a series for the line segments
    if get(d, :linewidth, 1) > 0
        @series begin
            xseg, yseg, zseg = Segments(), Segments(), Segments()
            for (si, di, wi) in zip(source, destiny, weights)
                # add a line segment
                xsi, ysi, xdi, ydi = shorten_segment(x[si], y[si], x[di], y[di], shorten)
                # ysi, ydi = shorten_segment(y[si], y[di])
                if curves
                    xpt, ypt = random_control_point(xsi, xdi,
                                                    ysi, ydi,
                                                    curvature_scalar)
                    push!(xseg, xsi, xpt, xdi)
                    push!(yseg, ysi, ypt, ydi)
                    _3d && push!(zseg, z[si], z[si], z[di])
                else
                    push!(xseg, xsi, xdi)
                    push!(yseg, ysi, ydi)
                    _3d && push!(zseg, z[si], z[di])
                end
            end

            # generate a list of colors, one per segment
            grad = get(d, :linecolor, nothing)
            if isa(grad, ColorGradient)
                line_z := weights
            end

            seriestype := (curves ? :curves : (_3d ? :path3d : :path))
            series_annotations := []
            linewidth --> 1
            markershape := :none
            markercolor := :black
            primary := false
            _3d ? (xseg.pts, yseg.pts, zseg.pts) : (xseg.pts, yseg.pts)
        end
    end

    xlims --> extrema_plus_buffer(x)
    ylims --> extrema_plus_buffer(y)
    if _3d
        zlims --> extrema_plus_buffer(z)
    end

    seriestype := (_3d ? :scatter3d : :scatter)
    linewidth := 0
    linealpha := 0
    foreground_color_border --> nothing
    grid --> false
    legend --> false
    ticks --> nothing
    series_annotations --> map(string,names)
    markersize --> 10 + 100node_weights / sum(node_weights)
    _3d ? (x, y, z) : (x, y)
end

function extrema_plus_buffer(v)
    vmin,vmax = extrema(v)
    vdiff = vmax-vmin
    buffer = vdiff * 0.1
    vmin - buffer, vmax + buffer
end

function shorten_segment(x1, y1, x2, y2, shorten)
    xshort = shorten * (x2-x1)
    yshort = shorten * (y2-y1)
    x1+xshort, y1+yshort, x2-xshort, y2-yshort
end

# ---------------------------------------------------------------------------
# Arc/Chord Diagrams


function arcvertices{T}(source::AVec{T}, destiny::AVec{T})
    values = unique(vcat(source, destiny))
    [(i, i) for i in values ]
end

function arcvertices{T<:Union{Char,Symbol,AbstractString}}(source::AVec{T}, 
                                                                   destiny::AVec{T})
    lab2x = Dict{T,Int}()
    n = 1
    for element in vcat(source, destiny)
        if !haskey(lab2x, element)
            lab2x[element] = n
            n += 1
        end
    end 
    lab2x
end

@userplot ArcDiagram

@recipe function f(h::ArcDiagram)
    
    source, destiny, weights = get_source_destiny_weight(h.args...)
    
    vertices = arcvertices(source, destiny)
    
    # Box setup
    legend --> false
    aspect_ratio --> :equal
    grid --> false
    foreground_color_axis --> nothing
    foreground_color_border --> nothing
    ticks --> nothing
    
    usegradient = length(unique(weights)) != 1

    if usegradient
        colorgradient = ColorGradient(get(d,:linecolor,cgrad()))
        wmin,wmax = extrema(weights)
    end

    for (i, j, value) in zip(source,destiny,weights)
        @series begin
            
            xi = vertices[i]
            xj = vertices[j]
            
            if usegradient
                linecolor --> colorgradient[(value-wmin)/(wmax-wmin)]
            end
            
            legend --> false
            label := ""
            primary := false

            r  = (xj - xi) / 2
            x₀ = (xi + xj) / 2
            θ = linspace(0,π,30)
            x₀ + r * cos(θ), r * sin(θ)
        end
    end
    
    @series begin
        if eltype(keys(vertices)) <: Union{Char, AbstractString, Symbol}
            series_annotations --> collect(keys(vertices))
            markersize --> 0
            y = -0.1
        else
            y =  0.0
        end
        seriestype --> :scatter
        collect(values(vertices)) , Float64[ y for i in vertices ]
    end
end


# =================================================
# Arc and chord diagrams

# ---------------------------------------------------------------------------
# Chord diagram

function arcshape(θ1, θ2)
    Plots.shape_coords(Shape(vcat(
        Plots.partialcircle(θ1, θ2, 15, 1.1),
        reverse(Plots.partialcircle(θ1, θ2, 15, 0.9))
    )))
end

# """
# `chorddiagram(source, destiny, weights[, grad, zcolor, group])`

# Plots a chord diagram, form `source` to `destiny`,
# using `weights` to determine the edge colors using `grad`.
# `zcolor` or `group` can be used to determine the node colors.
# """

@userplot ChordDiagram

@recipe function f(h::ChordDiagram)
    source, destiny, weights = get_source_destiny_weight(h.args...)

    xlims := (-1.2,1.2)
    ylims := (-1.2,1.2)
    legend := false
    grid := false
    xticks := nothing
    yticks := nothing

    nodemin, nodemax = extrema(vcat(source, destiny))
    weightmin, weightmax = extrema(weights)

    A  = 1.5π # Filled space
    B  = 0.5π # White space (empirical)

    Δα = A / nodemax
    Δβ = B / nodemax
    δ = Δα  + Δβ

    grad = cgrad(get(d,:linecolor,:inferno), get(d, :linealpha, nothing))

    for i in 1:length(source)
        # TODO: this could be a shape with varying width
        @series begin
            seriestype := :curves
            x := [cos((source[i ]-1)*δ + 0.5Δα), 0.0, cos((destiny[i]-1)*δ + 0.5Δα)]
            y := [sin((source[i ]-1)*δ + 0.5Δα), 0.0, sin((destiny[i]-1)*δ + 0.5Δα)]
            linecolor := grad[(weights[i] - weightmin) / (weightmax - weightmin)]
            primary := false
            ()
        end
    end

    for n in 0:(nodemax-1)
        sx, sy = arcshape(n*δ, n*δ + Δα)
        @series begin
            seriestype := :shape
            x := sx
            y := sy
            ()
        end
    end
end


