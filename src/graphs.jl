const _graph_funcs = Dict{Symbol,Any}(
    :spectral => spectral_graph,
    :stress => by_axis_local_stress_graph,
    :tree => tree_graph,
    :buchheim => buchheim_graph,
    :arcdiagram => arc_diagram,
    :chorddiagram => chord_diagram
)

const _graph_inputs = Dict{Symbol,Any}(
    :spectral => :adjmat,
    :stress => :adjmat,
    :tree => :sourcedestiny,
    :buchheim => :adjlist,
    :arcdiagram => :sourcedestiny,
    :chorddiagram => :sourcedestiny
)

function prepare_graph_inputs(method::Symbol, inputs...)
    input_type = get(_graph_inputs, method, :sourcedestiny)
    if input_type == :adjmat
        (get_adjacency_matrix(inputs...),)
    elseif input_type == :sourcedestiny
        get_source_destiny_weight(inputs...)
    elseif input_type == :adjlist
        (get_adjacency_list(inputs...),)
    end
end

# -----------------------------------------------------

function get_source_destiny_weight(mat::AbstractArray{T,2}) where T
    nrow, ncol = size(mat)    # rows are sources and columns are destinies

    nosymmetric = !issymmetric(mat) # plots only triu for symmetric matrices
    nosparse = !issparse(mat) # doesn't plot zeros from a sparse matrix

    L = length(mat)

    source  = Array{Int}(undef, L)
    destiny = Array{Int}(undef, L)
    weights  = Array{T}(undef, L)

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

function get_source_destiny_weight(source::AbstractVector, destiny::AbstractVector)
    if length(source) != length(destiny)
        throw(ArgumentError("Source and destiny must have the same length."))
    end
    source, destiny, Float64[ 1.0 for i in source ]
end

function get_source_destiny_weight(source::AbstractVector, destiny::AbstractVector, weights::AbstractVector)
    if !(length(source) == length(destiny) == length(weights))
        throw(ArgumentError("Source, destiny and weights must have the same length."))
    end
    source, destiny, weights
end

function get_source_destiny_weight(adjlist::AbstractVector{V}) where V<:AbstractVector{Int}
    source = Int[]
    destiny = Int[]
    for (i,l) in enumerate(adjlist)
        for j in l
            push!(source, i)
            push!(destiny, j)
        end
    end
    get_source_destiny_weight(source, destiny)
end

# -----------------------------------------------------

function get_adjacency_matrix(mat::AbstractMatrix)
    mat
end

function get_adjacency_matrix(source::AbstractVector{Int}, destiny::AbstractVector{Int}, weights::AbstractVector)
    n = max(maximum(source), maximum(destiny))
    Matrix(sparse(source, destiny, weights, n, n))
end

function get_adjacency_matrix(adjlist::AbstractVector{V}) where V<:AbstractVector{Int}
    s,d,w = get_source_destiny_weight(adjlist)
    get_adjacency_matrix(s, d, w)
end

# -----------------------------------------------------

function get_adjacency_list(mat::AbstractMatrix)
    get_adjacency_list(get_source_destiny_weight(mat))
end

function get_adjacency_list(source::AbstractVector{Int}, destiny::AbstractVector{Int}, weights::AbstractVector)
    n = max(maximum(source), maximum(destiny))
    adjlist = [Int[] for i=1:n]
    for (s,d) in zip(source,destiny)
        push!(adjlist[s], d)
    end
    adjlist
end

function get_adjacency_list(adjlist::AbstractVector{V}) where V<:AbstractVector{Int}
    adjlist
end

# -----------------------------------------------------

function make_symmetric(A::AbstractMatrix)
    A = copy(A)
    for i=1:size(A,1), j=(i+1):size(A,2)
        A[i,j] = A[j,i] = A[i,j]+A[j,i]
    end
    A
end

function compute_laplacian(adjmat::AbstractMatrix, node_weights::AbstractVector)
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

import LightGraphs

# TODO: so much wasteful conversion... do better
function estimate_distance(adjmat::AbstractMatrix)
    source, destiny, weights = get_source_destiny_weight(sparse(adjmat))

    g = LightGraphs.Graph(adjmat)
    dists = convert(Matrix{Float64}, hcat(map(i->LightGraphs.dijkstra_shortest_paths(g, i).dists, LightGraphs.vertices(g))...))
    tot = 0.0; cnt = 0
    for (i,d) in enumerate(dists)
        if d < 1e10
            tot += d
            cnt += 1
        end
    end
    avg = cnt > 0 ? tot / cnt : 1.0
    for (i,d) in enumerate(dists)
        if d > 1e10
            dists[i] = 3avg
        end
    end
    dists
end

function get_source_destiny_weight(g::LightGraphs.AbstractGraph)
    source = Vector{Int}()
    destiny =  Vector{Int}()
    sizehint!(source, LightGraphs.nv(g))
    sizehint!(destiny, LightGraphs.nv(g))
    for e in LightGraphs.edges(g)
      push!(source, LightGraphs.src(e))
      push!(destiny, LightGraphs.dst(e))
    end
    get_source_destiny_weight(source, destiny)
end

function get_adjacency_matrix(g::LightGraphs.Graph)
    get_adjacency_matrix(get_source_destiny_weight(g)...)
end

function get_adjacency_list(g::LightGraphs.AbstractGraph)
    g.fadjlist
end

# -----------------------------------------------------

# a graphplot takes in either an (N x N) adjacency matrix
#   note: you may want to pass node weights to markersize or marker_z
# A graph has N nodes where adj_mat[i,j] is the strength of edge i --> j.  (adj_mat[i,j]==0 implies no edge)

# NOTE: this is for undirected graphs... adjmat should be symmetric and non-negative

@userplot GraphPlot

@recipe function f(g::GraphPlot;
                   dim = 2,
                   free_dims = nothing,
                   T = Float64,
                   curves = true,
                   curvature_scalar = 0.2,
                   root = :top,
                   node_weights = nothing,
                   names = [],
                   fontsize = 7,
                   nodeshape = :hexagon,
                   nodesize = 1,
                   x = nothing,
                   y = nothing,
                   z = nothing,
                   method = :stress,
                   func = get(_graph_funcs, method, by_axis_local_stress_graph),
                   shorten = 0.0,
                   axis_buffer = 0.2,
                   layout_kw = Dict{Symbol,Any}()
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

    if isnothing(node_weights)
        node_weights = ones(n)
    end
    @assert length(node_weights) == n

    xyz = _3d ? (x,y,z) : (x,y)
    numnothing = count(isnothing, xyz)

    # do we want to compute coordinates?
    if numnothing > 0
        if isnothing(free_dims)
            # compute free_dims
            free_dims = findall(isnothing, xyz)
        end
        x, y, z = func(
            prepare_graph_inputs(method, source, destiny, weights)...;
            node_weights = node_weights,
            dim = dim,
            free_dims = free_dims,
            root = root,
            layout_kw...
        )
    end

    # reorient the points after root
    if root in (:left,:right)
        x,y = y,-x
    end
    if root == :left
        x,y = -x, y
    end
    if root == :bottom
        x,y = x, -y
    end

    # center and rescale to the widest of all dimensions
    xyz = _3d ? (x,y,z) : (x,y)

    if method == :arcdiagram
        xl, yl = arcdiagram_limits(x, source, destiny)
        xlims --> xl
        ylims --> yl
        ratio --> :equal
    elseif axis_buffer < 0 # equal axes
        ahw = 1.2 * 0.5 * maximum(v -> maximum(v)-minimum(v), xyz)
        xcenter = mean(extrema(x))
        #xlims --> (xcenter-ahw, xcenter+ahw)
        ycenter = mean(extrema(y))
        #ylims --> (ycenter-ahw, ycenter+ahw)
        if _3d
            zcenter = mean(extrema(z))
            #zlims --> (zcenter-ahw, zcenter+ahw)
        end
    else
        #xlims --> Plots.extrema_plus_buffer(x, axis_buffer)
        #ylims --> Plots.extrema_plus_buffer(y, axis_buffer)
        if _3d
            #zlims --> Plots.extrema_plus_buffer(z, axis_buffer)
        end
    end

    # create a series for the line segments
    if get(plotattributes, :linewidth, 1) > 0
        @series begin
            xseg = Vector{Float64}()
            yseg = Vector{Float64}()
            zseg = Vector{Float64}()
            l_wg = Vector{Float64}()
            for (si, di, wi) in zip(source, destiny, weights)

                # TO DO : Colouring edges by weight
                # add a line segment
                xsi, ysi, xdi, ydi = shorten_segment(x[si], y[si], x[di], y[di], shorten)
                if curves
                    if method in (:tree, :buchheim)
                        # for trees, shorten should be on one axis only
                        # dist = sqrt((x[di]-x[si])^2 + (y[di]-y[si])^2) * shorten
                        dist = shorten * (root in (:left,:bottom) ? 1 : -1)
                        ishoriz = root in (:left,:right)
                        xsi, xdi = (ishoriz ? (x[si]+dist,x[di]-dist) : (x[si],x[di]))
                        ysi, ydi = (ishoriz ? (y[si],y[di]) : (y[si]+dist,y[di]-dist))

                        xpts, ypts = directed_curve(xsi, xdi, ysi, ydi,
                                    xview=get(plotattributes, :xlims, (0,1)), yview=get(plotattributes, :ylims, (0,1)), root=root)

                        append!(xseg, push!(xpts, NaN))
                        append!(yseg, push!(ypts, NaN))
                        append!(l_wg, [ wi for i in 1:length(xpts) ] )
                    elseif method == :arcdiagram
                        r  = (xdi - xsi) / 2
                        x₀ = (xdi + xsi) / 2
                        θ = range(0, stop=π, length=30)
                        xpts = x₀ .+ r .* cos.(θ)
                        ypts = r .* sin.(θ) .+ ysi # ysi == ydi
                        for x in xpts
                            push!(xseg, x)
                            push!(l_wg, wi)
                        end
                        push!(xseg, NaN)
                        for y in ypts
                            push!(yseg, y)
                        end
                        push!(yseg, NaN)
                    else
                        xpt, ypt = if method != :chorddiagram
                            random_control_point(xsi, xdi,
                                                 ysi, ydi,
                                                 curvature_scalar)
                        else
                            (0.0, 0.0)
                        end
                        push!(xseg, xsi, xpt, xdi, NaN)
                        push!(yseg, ysi, ypt, ydi, NaN)
                        _3d && push!(zseg, z[si], z[si], z[di], NaN)
                        push!(l_wg, wi)
                    end
                else
                    push!(xseg, xsi, xdi, NaN)
                    push!(yseg, ysi, ydi, NaN)
                    _3d && push!(zseg, z[si], z[di], NaN)
                    push!(l_wg, wi)
                end
            end

            # generate a list of colors, one per segment
            grad = get(plotattributes, :linecolor, nothing)
            if isa(grad, ColorGradient)
                line_z := l_wg
            end

            seriestype := (curves ? :curves : (_3d ? :path3d : :path))
            series_annotations := nothing
            linewidth --> 1
            markershape := :none
            markercolor := :black
            primary := false
            _3d ? (xseg, yseg, zseg) : (xseg, yseg)
        end
    end

    framestyle := :none
    axis := nothing
    legend --> false

    if method == :chorddiagram
        seriestype := :scatter
        markersize := 0
        markeralpha := 0
        ratio --> :equal
        if length(names) == length(x)
            annotations := [(x[i], y[i], names[i]) for i in 1:length(x)]
        end
        @series begin
            seriestype := :shape
            N = length(x)
            angles = Vector{Float64}(undef, N)
            for i in 1:N
                if y[i] > 0
                    angles[i] = acos(x[i])
                else
                    angles[i] = 2pi - acos(x[i])
                end
            end
            δ = 0.4 * (angles[2] - angles[1])
            vec_vec_xy = [ arcshape(Θ-δ,Θ+δ) for Θ in angles ] # Shape
            [ [ xy[1] for xy in vec_xy ] for vec_xy in vec_vec_xy ], [ [ xy[2] for xy in vec_xy ] for vec_xy in vec_vec_xy ]
        end
    else
        if isempty(names)
            seriestype := (_3d ? :scatter3d : :scatter)
            linewidth := 0
            linealpha := 0
            series_annotations --> map(string,names)
            markersize --> (10 .+ (100 .* node_weights) ./ sum(node_weights))
        else
            @assert !_3d  # TODO: make this work in 3D
            scalefactor = pop!(plotattributes, :markersize, nodesize)
            seriestype := :scatter
            # markersize := nodesize
            nodeshape = get(plotattributes, :markershape, nodeshape)
            # nodeshape = if isa(nodeshape, AbstractArray)
            #    [(sym) for sym in nodeshape] # Shape
            # else
            #    nodeshape # Shape
            # end
            # font(fontsize)
            series_annotations := (map(string,names), nodeshape, fontsize, scalefactor)
        end
    end
    xyz
end
