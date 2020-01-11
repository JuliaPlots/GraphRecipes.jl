const _graph_funcs = Dict{Symbol,Any}(
    :spectral => spectral_graph,
    :sfdp => sfdp_graph,
    :circular => circular_graph,
    :shell => shell_graph,
    :spring => spring_graph,
    :stress => by_axis_local_stress_graph,
    :tree => tree_graph,
    :buchheim => buchheim_graph,
    :arcdiagram => arc_diagram,
    :chorddiagram => chord_diagram
)

const _graph_inputs = Dict{Symbol,Any}(
    :spectral => :adjmat,
    :sfdp => :adjmat,
    :circular => :adjmat,
    :shell => :adjmat,
    :stress => :adjmat,
    :spring => :adjmat,
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
    deg = vec(sum(adjmat; dims=1)) - diag(adjmat)
    D = diagm(0 => deg)

    # Laplacian (L = D - adjmat)
    L = eltype(adjmat)[i == j ? deg[i] : -adjmat[i,j] for i=1:n,j=1:n]

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

function get_adjacency_matrix(g::LightGraphs.AbstractGraph)
    adjacency_matrix(g)
end

function get_adjacency_matrix(source::AbstractVector{Int}, destiny::AbstractVector{Int})
    get_adjacency_matrix(source, destiny, ones(length(source)))
end

function get_adjacency_list(g::LightGraphs.AbstractGraph)
    g.fadjlist
end

# -----------------------------------------------------

# a graphplot takes in either an (N x N) adjacency matrix
#   note: you may want to pass node weights to markersize or marker_z
# A graph has N nodes where adj_mat[i,j] is the strength of edge i --> j.  (adj_mat[i,j]==0 implies no edge)

# NOTE: this is for undirected graphs... adjmat should be symmetric and non-negative

const graph_aliases = Dict(:curvature_scalar => [:curvaturescalar,:curvature],
                           :node_weights => [:nodeweights],
                           :nodeshape => [:node_shape],
                           :nodesize => [:node_size],
                           :nodecolor => [:marker_color],
                           :shorten => [:shorten_edge],
                           :axis_buffer => [:axisbuffer],
                           :edgewidth => [:edge_width,:ew],
                           :edgelabel => [:edge_label,:el],
                           :edgelabel_offset => [:edgelabeloffset,:elo],
                           :self_edge_size => [:selfedgesize,:ses],
                           :edge_label_box => [:edgelabelbox,:edgelabel_box,:elb],
)

@userplot GraphPlot

@recipe function f(g::GraphPlot;
                   dim = 2,
                   free_dims = nothing,
                   T = Float64,
                   curves = true,
                   curvature_scalar = 0.05,
                   root = :top,
                   node_weights = nothing,
                   names = [],
                   fontsize = 7,
                   nodeshape = :hexagon,
                   nodesize = 0.1,
                   nodecolor = 1,
                   x = nothing,
                   y = nothing,
                   z = nothing,
                   method = :stress,
                   func = get(_graph_funcs, method, by_axis_local_stress_graph),
                   shorten = 0.0,
                   axis_buffer = 0.2,
                   layout_kw = Dict{Symbol,Any}(),
                   edgewidth = (s,d,w)->1,
                   edgelabel = nothing,
                   edgelabel_offset = 0.0,
                   self_edge_size = 0.1,
                   edge_label_box = true,
                  )
    # Process GraphRecipes aliases, keywords that are already in Plots.jl have to be dealt
    # with individually, since they cannot be deleted from the plotattributes
    # dictionary.
    if haskey(plotattributes, :markercolor)
        plotattributes[:nodecolor] = nodecolor
        nodecolor = plotattributes[:markercolor]
        delete!(plotattributes, :nodecolor)
    end
    if haskey(plotattributes, :markersize)
        plotattributes[:nodesize] = nodesize
        nodesize = plotattributes[:markersize]
        delete!(plotattributes, :nodesize)
    end
    if haskey(plotattributes, :markershape)
        plotattributes[:nodeshape] = nodeshape
        nodeshape = plotattributes[:markershape]
        delete!(plotattributes, :nodeshape)
    end

    curvature_scalar = replace_kwarg!(curvature_scalar, :curvature_scalar, plotattributes, graph_aliases)
    node_weights = replace_kwarg!(node_weights, :node_weights, plotattributes, graph_aliases)
    nodeshape = replace_kwarg!(nodeshape, :nodeshape, plotattributes, graph_aliases)
    nodesize = replace_kwarg!(nodesize, :nodesize, plotattributes, graph_aliases)
    nodecolor = replace_kwarg!(nodecolor, :nodecolor, plotattributes, graph_aliases)
    shorten = replace_kwarg!(shorten, :shorten, plotattributes, graph_aliases)
    axis_buffer = replace_kwarg!(axis_buffer, :axis_buffer, plotattributes, graph_aliases)
    edgewidth = replace_kwarg!(edgewidth, :edgewidth, plotattributes, graph_aliases)
    edgelabel = replace_kwarg!(edgelabel, :edgelabel, plotattributes, graph_aliases)
    edgelabel_offset = replace_kwarg!(edgelabel_offset, :edgelabel_offset, plotattributes, graph_aliases)
    self_edge_size = replace_kwarg!(self_edge_size, :self_edge_size, plotattributes, graph_aliases)
    edge_label_box = replace_kwarg!(edge_label_box, :edge_label_box, plotattributes, graph_aliases)
    edgelabel = replace_kwarg!(edgelabel, :edgelabel, plotattributes, graph_aliases)

    @assert dim in (2, 3)
    _3d = dim == 3
    adj_mat = get_adjacency_matrix(g.args...)
    isdirected = (g.args[1] isa DiGraph || !issymmetric(adj_mat)) && !in(method, (:tree, :buchheim))
    if isdirected && (g.args[1] isa Matrix)
        g = GraphPlot((adjacency_matrix(DiGraph(g.args[1])),))
    end

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
        xlims = ignorenan_extrema(x)
        if method != :chorddiagram && numnothing > 0
            x .-= mean(x)
            x /= (xlims[2] - xlims[1])
            y .-= mean(y)
            ylims = ignorenan_extrema(y)
            y /= (ylims[2] - ylims[1])
        end
        xlims --> extrema_plus_buffer(x, axis_buffer)
        ylims --> extrema_plus_buffer(y, axis_buffer)
        if _3d
            if method != :chorddiagram && numnothing > 0
                zlims = ignorenan_extrema(z)
                z .-= mean(z)
                z /= (zlims[2] - zlims[1])
            end
            zlims --> extrema_plus_buffer(z, axis_buffer)
        end
    end
    # center and rescale to the widest of all dimensions
    xyz = _3d ? (x,y,z) : (x,y)
    # Get the coordinates for the edges of the nodes.
    node_vec_vec_xy = []
    nodewidth = 0.0
    nodewidth_array = Vector{Float64}(undef, length(x))
    if !_3d
        for i in 1:length(x)
            node_weight = isnothing(node_weights) ? 1 : (10 + 100node_weights[i]/sum(node_weights))/50
            xextent, yextent = if isempty(names)
                [x[i] .+ [-0.5nodesize*node_weight, 0.5nodesize*node_weight], y[i] .+ [-0.5nodesize*node_weight, 0.5nodesize*node_weight]]
            else
                annotation_extent(plotattributes,
                                  (x[i], y[i], names[ifelse(i % length(names) == 0, length(names),
                                                            i % length(names))], fontsize*nodesize*node_weight))
            end
            nodewidth = xextent[2] - xextent[1]
            nodewidth_array[i] = nodewidth
            if nodeshape == :circle
                push!(node_vec_vec_xy, partialcircle(0, 2π, [x[i], y[i]],
                                                     80, nodewidth/2))
            elseif (nodeshape == :rect) || (nodeshape == :rectangle)
                push!(node_vec_vec_xy, [(xextent[1],yextent[1]),
                                        (xextent[2],yextent[1]),
                                        (xextent[2],yextent[2]),
                                        (xextent[1],yextent[2]),
                                        (xextent[1],yextent[1])])
            elseif nodeshape == :hexagon
                push!(node_vec_vec_xy, partialcircle(0, 2π, [x[i], y[i]],
                                                     7, nodewidth/2))
            elseif nodeshape == :ellipse
                nodeheight = (yextent[2] - yextent[1])
                push!(node_vec_vec_xy, partialellipse(0, 2π, [x[i], y[i]],
                                                      80, nodewidth/2, nodeheight/2))
            else
                error("Unknown nodeshape: $(nodeshape). Choose from :circle, ellipse, :hexagon, :rect or :rectangle.")
            end
        end
    else
        @assert _3d # TODO Make 3d work.
    end
    # create a series for the line segments
    if get(plotattributes, :linewidth, 1) > 0
        # generate a list of colors, one per segment
        segment_colors = get(plotattributes, :linecolor, nothing)
        edge_label_array = Vector{Tuple}()
        edge_label_box_vertices_array = Vector{Array}()
        if !isa(edgelabel, Dict) && !isnothing(edgelabel)
            tmp = Dict()
            if length(size(edgelabel)) < 2
                matrix_size = round(Int, sqrt(length(edgelabel)))
                edgelabel = reshape(edgelabel, matrix_size, matrix_size)
            end
            for i in 1:size(edgelabel)[1]
                for j in 1:size(edgelabel)[2]
                    if islabel(edgelabel[i, j])
                        tmp[(i, j)] = edgelabel[i, j]
                    end
                end
            end
            edgelabel = tmp
        end
        # If the edgelabel dictionary is full of length two tuples, then make all of the
        # tuples length three with last element 1. (i.e. a multigraph that has no extra
        # edges).
        if edgelabel isa Dict
            edgelabel = convert(Dict{Any,Any}, edgelabel)
            for key in keys(edgelabel)
                if length(key) == 2
                    edgelabel[(key..., 1)] = edgelabel[key]
                end
            end
        end
        edge_has_been_seen = Dict()
        for edge in zip(source, destiny)
            edge_has_been_seen[edge] = 0
        end
        for (i, (si, di, wi)) in enumerate(zip(source, destiny, weights))
            edge_has_been_seen[(si, di)] += 1
            @series begin
                xseg = Vector{Float64}()
                yseg = Vector{Float64}()
                zseg = Vector{Float64}()
                l_wg = Vector{Float64}()

                # TO DO : Colouring edges by weight
                # add a line segment
                xsi, ysi, xdi, ydi = shorten_segment(x[si], y[si], x[di], y[di], shorten)
                # For directed graphs, shorten the line segment so that the edge ends at
                # the perimeter of the destiny node.
                θ = (edge_has_been_seen[(si, di)] - 1)*pi/8
                if isdirected && si != di && !_3d
                    α = atan(y[si] - y[di], x[si] - x[di])
                    if sign(si - di) < 0
                        α = x[si] < x[di] ? θ + α : α - θ
                    else
                        α = x[si] > x[di] ? θ + α : α - θ
                    end
                    xdi = x[di] + nodewidth_array[di]*cos(α)/2
                    ydi = y[di] + nodewidth_array[di]*sin(α)/2
                    if nodeshape != :circle
                        xsi, ysi, xdi, ydi = nearest_intersection(x[si], y[si], xdi, ydi,
                                                                  node_vec_vec_xy[di])
                    end
                end
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
                            control_point(xsi, x[di],
                                          ysi, y[di],
                                          edge_has_been_seen[(si, di)]*curvature_scalar*sign(si - di))
                        else
                            (0.0, 0.0)
                        end
                        xpts = [xsi, xpt, xdi]
                        ypts = [ysi, ypt, ydi]
                        t = range(0, stop=1, length=3)
                        A = hcat(xpts, ypts)

                        itp = scale(interpolate(A, BSpline(Cubic(Natural(OnGrid())))), t, 1:2)

                        tfine = range(0, stop=1, length=30)
                        xpts, ypts = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
                        if !isnothing(edgelabel) && haskey(edgelabel, (si, di, edge_has_been_seen[(si, di)]))
                            q = control_point(xsi, x[di],
                                              ysi, y[di],
                                              (edgelabel_offset
                                              + edge_has_been_seen[(si, di)]*curvature_scalar)*sign(si - di))
                            push!(edge_label_array,
                                  (q...,
                                   string(edgelabel[(si, di, edge_has_been_seen[(si, di)])]), fontsize))
                            edge_label_box_vertices = (
                            annotation_extent(plotattributes,
                                              (q[1], q[2],
                                              edgelabel[(si, di, edge_has_been_seen[(si, di)])],
                                              0.05fontsize)))
                            if !any(isnan.(q))
                                push!(edge_label_box_vertices_array, edge_label_box_vertices)
                            end
                        end
                        if method != :chorddiagram && !_3d
                            append!(xseg, push!(xpts, NaN))
                            append!(yseg, push!(ypts, NaN))
                            push!(l_wg, wi)
                        else
                            push!(xseg, xsi, xpt, xdi, NaN)
                            push!(yseg, ysi, ypt, ydi, NaN)
                            _3d && push!(zseg, z[si], z[si], z[di], NaN)
                            push!(l_wg, wi)
                        end
                    end
                else
                    push!(xseg, xsi, xdi, NaN)
                    push!(yseg, ysi, ydi, NaN)
                    _3d && push!(zseg, z[si], z[di], NaN)

                    if !isnothing(edgelabel) && haskey(edgelabel, (si, di, edge_has_been_seen[(si, di)]))
                        q = [(xsi + xdi)/2, (ysi + ydi)/2]
                        push!(edge_label_array,
                              (q...,
                               string(edgelabel[(si, di, edge_has_been_seen[(si, di)])]), fontsize))
                        edge_label_box_vertices = (
                        annotation_extent(plotattributes,
                                          (q[1], q[2],
                                          edgelabel[(si, di, edge_has_been_seen[(si, di)])],
                                          0.05fontsize)))
                        if !any(isnan.(q))
                            push!(edge_label_box_vertices_array, edge_label_box_vertices)
                        end
                    end
                end
            if si == di && !_3d
                inds = 1:n .!= si
                self_edge_angle = pi/8 + (edge_has_been_seen[(si, di)] - 1)*pi/8
                θ1 = unoccupied_angle(xsi, ysi, x[inds], y[inds]) - self_edge_angle/2
                θ2 = θ1 + self_edge_angle
                nodewidth = nodewidth_array[si]
                if nodeshape == :circle
                    xpts = [xsi + nodewidth*cos(θ1)/2,
                            NaN, NaN, NaN,
                            xsi + nodewidth*cos(θ2)/2]
                    xpts[2] = mean([xpts[1], xpts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*cos(θ1)
                    xpts[3] = mean([xpts[1], xpts[end]]) + edge_has_been_seen[(si, di)]*self_edge_size*cos((θ1 + θ2)/2)
                    xpts[4] = mean([xpts[1], xpts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*cos(θ2)
                    ypts = [ysi + nodewidth*sin(θ1)/2,
                            NaN, NaN, NaN,
                            ysi + nodewidth*sin(θ2)/2]
                    ypts[2] = mean([ypts[1], ypts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*sin(θ1)
                    ypts[3] = mean([ypts[1], ypts[end]]) + edge_has_been_seen[(si, di)]*self_edge_size*sin((θ1 + θ2)/2)
                    ypts[4] = mean([ypts[1], ypts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*sin(θ2)

                    t = range(0, stop=1, length=5)
                    A = hcat(xpts, ypts)

                    itp = scale(interpolate(A, BSpline(Cubic(Natural(OnGrid())))), t, 1:2)

                    tfine = range(0, stop=1, length=50)
                    xpts, ypts = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
                else
                    _, _,
                    start_point1,
                    start_point2 = nearest_intersection(xsi, ysi,
                                                        xsi + 2nodewidth*cos(θ1),
                                                        ysi + 2nodewidth*sin(θ1),
                                                        node_vec_vec_xy[si])
                    _, _, end_point1,
                    end_point2 = nearest_intersection(xsi + edge_has_been_seen[(si, di)]*(nodewidth + self_edge_size)*cos(θ2),
                                                      ysi + edge_has_been_seen[(si, di)]*(nodewidth + self_edge_size)*sin(θ2),
                                                      xsi,
                                                      ysi,
                                                      node_vec_vec_xy[si])

                    xpts = [start_point1,
                            NaN, NaN, NaN,
                            end_point1]
                    xpts[2] = mean([xpts[1], xpts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*cos(θ1)
                    xpts[3] = mean([xpts[1], xpts[end]]) + edge_has_been_seen[(si, di)]*self_edge_size*cos((θ1 + θ2)/2)
                    xpts[4] = mean([xpts[1], xpts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*cos(θ2)
                    ypts = [start_point2,
                            NaN, NaN, NaN,
                            end_point2]
                    ypts[2] = mean([ypts[1], ypts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*sin(θ1)
                    ypts[3] = mean([ypts[1], ypts[end]]) + edge_has_been_seen[(si, di)]*self_edge_size*sin((θ1 + θ2)/2)
                    ypts[4] = mean([ypts[1], ypts[end]]) + 0.5*(0.5 + edge_has_been_seen[(si, di)])*self_edge_size*sin(θ2)

                    t = range(0, stop=1, length=5)
                    A = hcat(xpts, ypts)

                    itp = scale(interpolate(A, BSpline(Cubic(Natural(OnGrid())))), t, 1:2)

                    tfine = range(0, stop=1, length=50)
                    xpts, ypts = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
                end
                append!(xseg, push!(xpts, NaN))
                append!(yseg, push!(ypts, NaN))
                mid_ind = div(length(xpts), 2)
                q = [xpts[mid_ind] + edgelabel_offset*cos((θ1 + θ2)/2),
                     ypts[mid_ind] + edgelabel_offset*sin((θ1 + θ2)/2)]
                if !isnothing(edgelabel) && haskey(edgelabel, (si, di, edge_has_been_seen[(si, di)]))
                    push!(edge_label_array,
                          (q...,
                           string(edgelabel[(si, di, edge_has_been_seen[(si, di)])]), fontsize))
                    edge_label_box_vertices = annotation_extent(plotattributes, (q...,
                                                                edgelabel[(si, di, edge_has_been_seen[(si, di)])],
                                                                0.05fontsize))
                    if !any(isnan.(q))
                        push!(edge_label_box_vertices_array, edge_label_box_vertices)
                    end
                end
            end

            if isa(segment_colors, ColorGradient)
                line_z := segment_colors[i]
            end
            linewidthattr = get(plotattributes, :linewidth, 1)
            seriestype := if method in (:tree, :buchheim, :chorddiagram)
                :curves
            else
                if _3d
                    if curves
                        :curves
                    else
                        :path3d
                    end
                else
                    :path
                end
            end
            linewidth --> linewidthattr * edgewidth(si, di, wi)
            markershape := :none
            markercolor := :black
            isdirected && (arrow --> :simple, :head, 0.3, 0.3)
            primary := false
            _3d ? (xseg, yseg, zseg) : (xseg, yseg)
            end

        end
    end
    # The boxes around edge labels are defined as another list of series that sits on top
    # of the series for the edges.
    edge_has_been_seen = Dict()
    for edge in zip(source, destiny)
        edge_has_been_seen[edge] = 0
    end
    if edge_label_box
        index = 0
        for (i, (si, di, wi)) in enumerate(zip(source, destiny, weights))
            edge_has_been_seen[(si, di)] += 1
            if !isnothing(edgelabel) && haskey(edgelabel, (si, di, edge_has_been_seen[(si, di)]))
                index += 1
                @series begin
                    seriestype := :shape
                    fillcolor --> get(plotattributes, :background_color, :white)
                    linewidth := 0
                    linealpha := 0
                    edge_label_box_vertices = edge_label_box_vertices_array[index]
                    ([edge_label_box_vertices[1][1], edge_label_box_vertices[1][2],
                      edge_label_box_vertices[1][2], edge_label_box_vertices[1][1],
                      edge_label_box_vertices[1][1]],
                     [edge_label_box_vertices[2][1], edge_label_box_vertices[2][1],
                      edge_label_box_vertices[2][2], edge_label_box_vertices[2][2],
                      edge_label_box_vertices[2][1]])
                end
            end
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
            if _3d
                seriestype := :scatter3d
                linewidth := 0
                linealpha := 0
                series_annotations --> map(string,names)
                markersize --> (10 .+ (100 .* node_weights) ./ sum(node_weights))
            else
                for (i, vec_xy) in enumerate(node_vec_vec_xy)
                    @series begin
                        if !_3d
                            seriestype := :shape
                            tmp = nodecolor isa AbstractArray ? nodecolor[i] : nodecolor
                            fillcolor --> tmp
                            linewidth := 0
                            linealpha := 0
                            ([xy[1] for xy in vec_xy], [xy[2] for xy in vec_xy])
                        else  # TODO make 3d work.
                            seriestype := :volume
                            ([[xyz[1] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz],
                             [[xyz[2] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz],
                             [[xyz[3] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz])
                        end
                    end
                    seriestype := :scatter
                    markersize := 0
                    markeralpha := 0
                end
                !isnothing(edgelabel) && (annotation := edge_label_array)
            end
        else
            @assert !_3d  # TODO: make this work in 3D
            for (i, vec_xy) in enumerate(node_vec_vec_xy)
                @series begin
                    if !_3d
                        seriestype := :shape
                        tmp = nodecolor isa AbstractArray ? nodecolor[i] : nodecolor
                        fillcolor --> tmp
                        linewidth := 0
                        linealpha := 0
                        ([xy[1] for xy in vec_xy], [xy[2] for xy in vec_xy])
                    else  # TODO make 3d work.
                        seriestype := :volume
                        ([[xyz[1] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz],
                         [[xyz[2] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz],
                         [[xyz[3] for xyz in vec_xyz] for vec_xyz in node_vec_vec_xyz])
                    end
                end
            end
            seriestype := :scatter
            markersize := 0
            markeralpha := 0
            annotations := [edge_label_array ; [(x[i], y[i],
                                                names[ifelse(i % length(names) == 0,
                                                              length(names),
                                                              i % length(names))],
                                                fontsize) for i in 1:length(x)]]
        end
    end
    xyz
end

@recipe function f(g::AbstractGraph)
    GraphPlot(get_source_destiny_weight(get_adjacency_list(g)))
end
