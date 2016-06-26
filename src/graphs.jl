
# a graphplot takes in either an (N x N) adjacency matrix
#   note: you may want to pass node weights to markersize or marker_z
# A graph has N nodes where adj_mat[i,j] is the strength of edge i --> j.  (adj_mat[i,j]==0 implies no edge)

# NOTE: this is for undirected graphs... adjmat should be symmetric and non-negative

@userplot GraphPlot

# see: http://www.research.att.com/export/sites/att_labs/groups/infovis/res/legacy_papers/DBLP-journals-camwa-Koren05.pdf
# also: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.3.2055&rep=rep1&type=pdf

# this recipe uses the technique of Spectral Graph Drawing, which is an
# under-appreciated method of graph layouts; easier, simpler, and faster
# than the more common spring-based methods.
function spectral_graph(adjmat::AbstractMatrix, node_weights::AbstractVector = ones(size(adjmat,1)))
    n, m = size(adjmat)
    @assert n == m == length(node_weights)

    # scale the edge values by the product of node_weights, so that "heavier" nodes also form
    # stronger connections
    adjmat = adjmat .* sqrt(node_weights * node_weights')

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


# -----------------------------------------------------

# TODO: maybe also implement Catmull-Rom Splines? http://www.mvps.org/directx/articles/catmull/

# -----------------------------------------------------

# # get the value of the curve point at position t 
# function bezier_value(pts::AVec, t::Real)
#     val = 0.0
#     n = length(pts)-1
#     for (i,p) in enumerate(pts)
#         val += p * binomial(n, i-1) * (1-t)^(n-i+1) * t^(i-1)
#     end
#     val
# end

# # create segmented bezier curves in place of line segments
# @recipe function f(::Type{Val{:curves}}, x, y, z)
#     _3d = z != nothing
#     args = _3d ? (x,y,z) : (x,y)
#     newx, newy = zeros(0), zeros(0)
#     newz = _3d ? zeros(0) : nothing
#     for rng in iter_segments(args...)
#         length(rng) < 2 && continue
#         ts = linspace(0, 1, pop!(d, :npoints, 20))
#         Plots.nanappend!(newx, map(t -> bezier_value(x[rng],t), ts))
#         Plots.nanappend!(newy, map(t -> bezier_value(y[rng],t), ts))
#         if _3d
#             Plots.nanappend!(newz, map(t -> bezier_value(z[rng],t), ts))
#         end
#         @show rng ts newx newy newz
#     end
#     seriestype := :path
#     x := newx
#     y := newy
#     z := newz
#     ()
# end

# # convert pairs of ((x1,y1), (x2,y2)) or ((x1,y1,z2), (x2,y2,z2)) into valid line segments
# @recipe function f(::Type{Val{:segments}}, x, y, z)
#     # for each segment, we need 3 values in each of x/y/z/line_z
#     _3d = z != nothing
#     n = max(length(x), length(y))
#     lx, ly = zeros(3n), zeros(3n)
#     lz = _3d ? zeros(3n) : nothing
#     # lx, ly, lz = zeros(T,0), zeros(T,0), zeros(T,0)
#     # line_z = zeros(T,0)
#     line_z = zeros(3n)
#     for
#     for i=1:n, j=i+1:n
#         aij = adjmat[i,j]
#         if aij ≉ 0
#             # @show aij, x,y,i,j
#             append!(lx, T[x[i], x[j], NaN])
#             append!(ly, T[y[i], y[j], NaN])
#             if dim == 3
#                 append!(lz, T[z[i], z[j], NaN])
#             end
#             append!(line_z, T[aij, aij])
#             # TODO: when supported, add line width for this line segment
#         end
#     end
# end

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


@recipe function f(g::GraphPlot; dim = 2,
                                 T = Float64,
                                 curves = true,
                                 curvature_scalar = 1,
                                 func = spectral_graph)
    @assert dim in (2, 3)
    delete!(d, :dim)
    delete!(d, :T)
    delete!(d, :curves)
    delete!(d, :curvature_scalar)
    delete!(d, :func)
    _3d = dim == 3

    adjmat = g.args[1]
    n, m = size(adjmat)
    x, y, z = func(g.args...)

    # create a series for the line segments
    if get(d, :linewidth, 1) > 0
        @series begin
            lx, ly = zeros(0), zeros(0)
            lz = _3d ? zeros(0) : nothing
            line_z = zeros(0)

            # skipped when user overrides linewidth to 0
            # we want to build new lx/ly/lz for the lines
            # note: we only do the lower triangle
            for i=2:n, j=1:i-1
                aij = adjmat[i,j]
                if aij ≉ 0
                    # @show aij, x,y,i,j
                    # add the first point, NaN-separated
                    Plots.nanpush!(lx, x[i])
                    Plots.nanpush!(ly, y[i])
                    _3d && Plots.nanpush!(lz, z[i])

                    # add curve control points?
                    if curves
                        xpt, ypt = random_control_point(x[i], x[j], y[i], y[j], curvature_scalar)
                        push!(lx, xpt)
                        push!(ly, ypt)
                        # @show (x[i], xpt, x[j]), (y[i], ypt, y[j])
                        _3d && push!(lz, 0.5(z[i] + z[j]))
                    end

                    # add the last point and line_z value
                    push!(lx, x[j])
                    push!(ly, y[j])
                    _3d && push!(lz, z[j])
                    push!(line_z, aij)
                end
            end

            # update line_z to the correct size
            if isa(get(d, :linecolor, nothing), ColorGradient)
                line_z = vec(repmat(line_z', curves ? 4 : 3, 1))
                line_z --> line_z, :quiet
            end

            seriestype := (curves ? :curves : (_3d ? :path : :path3d))
            series_annotations := []
            # linecolor --> :black
            linewidth --> 1
            markershape := :none
            markercolor := :black
            primary := false
            # Plots.DD(d)
            _3d ? (lx, ly, lz) : (lx, ly)
        end
    end

    seriestype := (_3d ? :scatter3d : :scatter)
    linewidth := 0
    linealpha := 0
    foreground_color_border --> nothing
    grid --> false
    legend --> false
    ticks --> nothing
    if length(g.args) > 1
        node_weights = g.args[2]
        markersize --> 10 + 100node_weights / sum(node_weights)
    end
    _3d ? (x, y, z) : (x, y)
end

# ---------------------------------------------------------------------------
# Arc/Chord Diagrams

function get_source_destiny_weight{T}(mat::AbstractArray{T,2})
    nrow, ncol = size(mat)    # rows are sources and columns are destinies

    nosymmetric = !issym(mat) # plots only triu for symmetric matrices
    nosparse = !issparse(mat) # doesn't plot zeros from a sparse matrix

    L = length(mat)

    source  = Array(Int, L)
    destiny = Array(Int, L)
    weight  = Array(T, L)

    idx = 1
    for i in 1:nrow, j in 1:ncol
        value = mat[i, j]
        if !isnan(value) && ( nosparse || value != zero(T) ) # TODO: deal with Nullable

            if i < j
                source[idx]  = i
                destiny[idx] = j
                weight[idx]  = value
                idx += 1
            elseif nosymmetric && (i > j)
                source[idx]  = i
                destiny[idx] = j
                weight[idx]  = value
                idx += 1
            end

        end
    end

    resize!(source, idx-1), resize!(destiny, idx-1), resize!(weight, idx-1)
end

function get_source_destiny_weight(source::AbstractVector, destiny::AbstractVector)
    if length(source) != length(destiny)
        throw(ArgumentError("Source and destiny must have the same length."))
    end
    source, destiny, Float64[ 1.0 for i in source ]
end

function get_source_destiny_weight(source::AbstractVector, destiny::AbstractVector, weight::AbstractVector)
    if !(length(source) == length(destiny) == length(weight))
        throw(ArgumentError("Source, destiny and weight must have the same length."))
    end
    source, destiny, weight
end

function arcvertices{T}(source::AbstractVector{T}, destiny::AbstractVector{T})
    values = unique(vcat(source, destiny))
    [ i => i for i in values ]
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
    
    source, destiny, weight = get_source_destiny_weight(h.args...)
    
    vertices = arcvertices(source, destiny)
    
    # Box setup
    legend --> false
    aspect_ratio --> :equal
    grid --> false
    foreground_color_axis --> nothing
    foreground_color_border --> nothing
    ticks --> nothing
    
    usegradient = length(unique(weight)) != 1

    if usegradient
        colorgradient = ColorGradient(get(d,:linecolor,Plots.default_gradient()))
        wmin,wmax = extrema(weight)
    end

    for (i, j, value) in zip(source,destiny,weight)
        @series begin
            
            xi = vertices[i]
            xj = vertices[j]
            
            if usegradient
                linecolor --> getColorZ(colorgradient, (value-wmin)/(wmax-wmin))
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

# "Takes an adjacency matrix and returns source, destiny and weight lists"
# function mat2list{T}(mat::AbstractArray{T,2})
#     nrow, ncol = size(mat) # rows are sources and columns are destinies

#     nosymmetric = !issym(mat) # plots only triu for symmetric matrices
#     nosparse = !issparse(mat) # doesn't plot zeros from a sparse matrix

#     L = length(mat)

#     source  = Array(Int, L)
#     destiny = Array(Int, L)
#     weight  = Array(T, L)

#     idx = 1
#     for i in 1:nrow, j in 1:ncol
#         value = mat[i, j]
#         if !isnan(value) && ( nosparse || value != zero(T) ) # TODO: deal with Nullable

#             if i < j
#                 source[idx]  = i
#                 destiny[idx] = j
#                 weight[idx]  = value
#                 idx += 1
#             elseif nosymmetric && (i > j)
#                 source[idx]  = i
#                 destiny[idx] = j
#                 weight[idx]  = value
#                 idx += 1
#             end

#         end
#     end

#     resize!(source, idx-1), resize!(destiny, idx-1), resize!(weight, idx-1)
# end

# curvecolor(value, min, max, grad) = getColorZ(grad, (value-min)/(max-min))

# ---------------------------------------------------------------------------
# Chord diagram

function arcshape(θ1, θ2)
    Plots.shape_coords(Shape(vcat(
        Plots.partialcircle(θ1, θ2, 15, 1.1),
        reverse(Plots.partialcircle(θ1, θ2, 15, 0.9))
    )))
end

# colorlist(grad, ::Void) = :darkgray

# function colorlist(grad, z)
#     zmin, zmax = extrema(z)
#     RGBA{Float64}[getColorZ(grad, (zi-zmin)/(zmax-zmin)) for zi in z]'
# end

# """
# `chorddiagram(source, destiny, weight[, grad, zcolor, group])`

# Plots a chord diagram, form `source` to `destiny`,
# using `weight` to determine the edge colors using `grad`.
# `zcolor` or `group` can be used to determine the node colors.
# """
# function chorddiagram(source, destiny, weight; kargs...)


@userplot ChordDiagram

@recipe function f(h::ChordDiagram)
    
    source, destiny, weight = get_source_destiny_weight(h.args...)

    # args  = KW(kargs)
    # grad  = pop!(args, :grad,   ColorGradient([colorant"darkred", colorant"darkblue"]))
    # zcolor= pop!(args, :zcolor, nothing)
    # group = pop!(args, :group,  nothing)

    # if zcolor !== nothing && group !== nothing
    #     throw(ErrorException("group and zcolor can not be used together."))
    # end

    # if length(source) == length(destiny) == length(weight)
    xlims := (-1.2,1.2)
    ylims := (-1.2,1.2)
    legend := false
    grid := false
    xticks := nothing
    yticks := nothing

    nodemin, nodemax = extrema(vcat(source, destiny))
    weightmin, weightmax = extrema(weight)

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
            linecolor := grad[(weight[i] - weightmin) / (weightmax - weightmin)]
            primary := false
            ()
        end
        # curve = BezierCurve(P2[ (cos((source[i ]-1)*δ + 0.5Δα), sin((source[i ]-1)*δ + 0.5Δα)), (0,0),
        #                         (cos((destiny[i]-1)*δ + 0.5Δα), sin((destiny[i]-1)*δ + 0.5Δα)) ])
        # plot!(curve_points(curve), line = grad[(weight[i] - weightmin) / (weightmax - weightmin)], 1, 1))
    end

    # if group === nothing
    #     c =  colorlist(grad, zcolor)
    # elseif length(group) == nodemax

    #     idx = collect(0:(nodemax-1))

    #     for g in group
    #         plot!([arcshape(n*δ, n*δ + Δα) for n in idx[group .== g]]; args...)
    #     end

    #     return plt

    # else
    #     throw(ErrorException("group should the ", nodemax, " elements."))
    # end

    # grad = cgrad(d[:fillcolor], d[:fillalpha])

    for n in 0:(nodemax-1)
        sx, sy = arcshape(n*δ, n*δ + Δα)
        @series begin
            seriestype := :shape
            x := sx
            y := sy
            ()
        end
    end

    # plot!([arcshape(n*δ, n*δ + Δα) for n in 0:(nodemax-1)]; mc=c, args...)

    # return plt

    # else
    #     throw(ArgumentError("source, destiny and weight should have the same length"))
    # end
end

# """
# `chorddiagram(mat[, grad, zcolor, group])`

# Plots a chord diagram from an adjacency matrix,
# using the values on the matrix as weights to determine edge colors.
# Doesn't show edges with value zero if the input is sparse.
# For simmetric matrices, only the upper triangular values are used.
# `zcolor` or `group` can be used to determine the node colors.
# """
# chorddiagram(mat::AbstractMatrix; kargs...) = chorddiagram(mat2list(mat)...; kargs...)

