
# -----------------------------------------------------
# -----------------------------------------------------

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

function spectral_graph(adjmat::AbstractMatrix; node_weights::AbstractVector = ones(size(adjmat,1)), kw...)
    positions = NetworkLayout.Spectral.layout(adjmat; node_weights=convert(Vector{Float64}, node_weights))

    ([p[1] for p in positions], [p[2] for p in positions], [p[3] for p in positions])
end

function spectral_graph(source::AbstractVector{Int}, destiny::AbstractVector{Int}, weights::AbstractVector; kw...)
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
# function by_axis_stress_graph(adjmat::AbstractMatrix, node_weights::AbstractVector = ones(size(adjmat,1));
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
function by_axis_local_stress_graph(adjmat::AbstractMatrix;
                              node_weights::AbstractVector = ones(size(adjmat,1)),
                              dim = 2, free_dims = 1:dim,
                              x = rand(length(node_weights)),
                              y = rand(length(node_weights)),
                              z = rand(length(node_weights)),
                              maxiter = 1000,
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
        if abs(thisstress - laststress) / abs(laststress) < 1e-6
            # info("converged. numiter=$k last=$laststress this=$thisstress")
            break
        end
        laststress = thisstress
    end

    dim == 2 ? (X..., nothing) : X
end

function by_axis_local_stress_graph(source::AbstractVector{Int}, destiny::AbstractVector{Int}, weights::AbstractVector; kw...)
    by_axis_local_stress_graph(get_adjacency_matrix(source, destiny, weights); kw...)
end

# -----------------------------------------------------

function buchheim_graph(adjlist::AbstractVector;
                    node_weights::AbstractVector = ones(length(adjlist)),
                    root::Symbol = :top,  # flow of tree: left, right, top, bottom
                    layers_scalar = 1.0,
                    layers = nothing,
                    dim = 2,
                    kw...)
    # @show adjlist typeof(adjlist)
    positions = NetworkLayout.Buchheim.layout(adjlist, nodesize = convert(Vector{Float64}, node_weights))
    Float64[p[1] for p in positions], Float64[p[2] for p in positions], nothing
end

# -----------------------------------------------------

function tree_graph(adjmat::AbstractMatrix; kw...)
    tree_graph(get_source_destiny_weight(adjmat)...; kw...)
end

function tree_graph(source::AbstractVector{Int}, destiny::AbstractVector{Int}, weights::AbstractVector;
                    node_weights::AbstractVector = ones(max(maximum(source), maximum(destiny))),
                    root::Symbol = :top,  # flow of tree: left, right, top, bottom
                    layers_scalar = 1.0,
                    layers = nothing,
                    positions = nothing,
                    dim = 2,
                    add_noise = true,
                    kw...)
    extrakw = Dict{Symbol,Any}(kw)
    # @show root layers positions dim add_noise extrakw
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
    if add_noise
        layers = layers + 0.6rand(size(layers)...)
    end

    # TODO: normalize layers somehow so it's in line with distances
    layers .*= layers_scalar
    if dim == 2
        if root in (:top, :bottom)
            extrakw[:y] = layers
            extrakw[:free_dims] = if isnothing(positions)
                [1]
            else
                extrakw[:x] = positions
                Int[]
            end
        elseif root in (:left, :right)
            extrakw[:x] = layers
            # extrakw[:free_dims] = [2]
            extrakw[:free_dims] = if isnothing(positions)
                [2]
            else
                extrakw[:y] = positions
                Int[]
            end
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

function arc_diagram(source::AbstractVector{Int},
                     destiny::AbstractVector{Int},
                     weights::AbstractVector;
                     kw...)
    X = unique(vcat(source, destiny))
    O = zeros(Int,length(X))
    X, O, O
end

# -----------------------------------------------------

function chord_diagram( source::AbstractVector{Int},
                        destiny::AbstractVector{Int},
                        weights::AbstractVector;
                        kw... )
    nodes = unique(vcat(source, destiny))
    N = length(nodes)
    δ = 2pi / N

    x = Array{Float64}(undef, N)
    y = Array{Float64}(undef, N)
    for i in 1:N
        x[i] = sin((i-1)*δ)
        y[i] = cos((i-1)*δ)
    end

    x, y, zeros(Int,length(x))
end
