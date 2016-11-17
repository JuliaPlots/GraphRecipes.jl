
include("graph_layouts.jl")

const _graph_funcs = KW(
    :spectral => spectral_graph,
    :stress => by_axis_local_stress_graph,
    :tree => tree_graph,
    :buchheim => buchheim_graph,
)

const _graph_inputs = KW(
    :spectral => :adjmat,
    :stress => :adjmat,
    :tree => :sourcedestiny,
    :buchheim => :adjlist,
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

function get_source_destiny_weight{V<:AbstractVector{Int}}(adjlist::AbstractVector{V})
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
    full(sparse(source, destiny, weights, n, n))
end

function get_adjacency_matrix{V<:AbstractVector{Int}}(adjlist::AbstractVector{V})
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

function get_adjacency_list{V<:AbstractVector{Int}}(adjlist::AbstractVector{V})
    adjlist
end

# -----------------------------------------------------

function make_symmetric(A::AbstractMatrix)
    A = copy(A)
    for i=1:size(A,1), j=i+1:size(A,2)
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

if Plots.is_installed("LightGraphs")
    @eval begin
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

        function get_source_destiny_weight(g::LightGraphs.SimpleGraph)
            source = Vector{Int}()
            destiny =  Vector{Int}()
            sizehint!(source, LightGraphs.nv(g))
            sizehint!(destiny, LightGraphs.nv(g))
            for (u, v) in LightGraphs.edges(g)
              push!(source, u)
              push!(destiny, v)
            end
            get_source_destiny_weight(source, destiny)
        end

        function get_adjacency_matrix(g::LightGraphs.SimpleGraph)
            get_adjacency_matrix(get_source_destiny_weight(g)...)
        end

        function get_adjacency_list(g::LightGraphs.SimpleGraph)
            g.fadjlist
        end
    end
else
    @eval function estimate_distance(adjmat::AbstractMatrix)
        warn("Install LightGraphs for the best layout calculations.")
        map(a -> a==0 ? 5.0 : a, adjmat)
    end
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
                   layout_kw = KW()
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
            free_dims = find(isnothing, xyz)
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

    # center and rescale to the widest of all dimensions
    xyz = _3d ? (x,y,z) : (x,y)

    if axis_buffer < 0 # equal axes
        ahw = 1.2 * 0.5 * maximum(v -> maximum(v)-minimum(v), xyz)
        xcenter = mean(extrema(x))
        xlims --> (xcenter-ahw, xcenter+ahw)
        ycenter = mean(extrema(y))
        ylims --> (ycenter-ahw, ycenter+ahw)
        if _3d
            zcenter = mean(extrema(z))
            zlims --> (zcenter-ahw, zcenter+ahw)
        end
    else
        xlims --> Plots.extrema_plus_buffer(x, axis_buffer)
        ylims --> Plots.extrema_plus_buffer(y, axis_buffer)
        if _3d
            zlims --> Plots.extrema_plus_buffer(z, axis_buffer)
        end
    end

    # create a series for the line segments
    if get(d, :linewidth, 1) > 0
        @series begin
            xseg, yseg, zseg = Segments(), Segments(), Segments()
            for (si, di, wi) in zip(source, destiny, weights)
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
                                    xview=d[:xlims], yview=d[:ylims], root=root)
                        push!(xseg, xpts)
                        push!(yseg, ypts)
                    else
                        xpt, ypt = random_control_point(xsi, xdi,
                                                        ysi, ydi,
                                                        curvature_scalar)
                        push!(xseg, xsi, xpt, xdi)
                        push!(yseg, ysi, ypt, ydi)
                        _3d && push!(zseg, z[si], z[si], z[di])
                    end
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
            series_annotations := nothing
            linewidth --> 1
            markershape := :none
            markercolor := :black
            primary := false
            _3d ? (xseg.pts, yseg.pts, zseg.pts) : (xseg.pts, yseg.pts)
        end
    end

    axis := nothing
    legend --> false

    if isempty(names)
        seriestype := (_3d ? :scatter3d : :scatter)
        linewidth := 0
        linealpha := 0
        series_annotations --> map(string,names)
        markersize --> 10 + 100node_weights / sum(node_weights)
    else
        @assert !_3d  # TODO: make this work in 3D
        seriestype := :scatter
        scalefactor = pop!(d, :markersize, nodesize)
        # markersize := nodesize
        nodeshape = get(d, :markershape, nodeshape)
        nodeshape = if isa(nodeshape, AbstractArray)
            [Shape(sym) for sym in nodeshape]
        else
            Shape(nodeshape)
        end
        series_annotations := (map(string,names), nodeshape, font(fontsize), scalefactor)
    end
    xyz
end

# ---------------------------------------------------------------------------
# Arc/Chord Diagrams


function arcvertices{T}(source::AbstractVector{T}, destiny::AbstractVector{T})
    values = unique(vcat(source, destiny))
    [(i, i) for i in values ]
end

function arcvertices{T<:Union{Char,Symbol,AbstractString}}(source::AbstractVector{T},
                                                                   destiny::AbstractVector{T})
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
    coords(Shape(vcat(
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
