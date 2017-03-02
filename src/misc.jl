
# ----------------------------------------------------------------
# Shapefile.jl plotting

# Remove legend, axes, grid, and force the aspect ratio for shape plots
@recipe function f(::Type{Val{:shapefile}}, x, y, z)
    legend --> false
    ticks --> nothing
    grid --> false
    aspect_ratio --> 1
    seriestype := :shape
    ()
end

if Plots.is_installed("Shapefile")
    @eval begin
        import Shapefile
        function shapefile_coords(poly::Shapefile.ESRIShape)
            start_indices = poly.parts+1
            end_indices = vcat(poly.parts[2:end], length(poly.points))
            x, y = zeros(0), zeros(0)
            for (si,ei) in zip(start_indices, end_indices)
                push!(x, NaN)
                push!(y, NaN)
                for pt in poly.points[si:ei]
                    push!(x, pt.x)
                    push!(y, pt.y)
                end
            end
            x, y
        end

        function shapefile_coords{T<:Shapefile.ESRIShape}(polys::AbstractArray{T})
            x, y = [], []
            for poly in polys
                xpart, ypart = shapefile_coords(poly)
                push!(x, xpart)
                push!(y, ypart)
            end
            x, y
        end

        @recipe f(poly::Shapefile.ESRIShape) = (seriestype --> :shapefile; shapefile_coords(poly))
        @recipe f{T<:Shapefile.ESRIShape}(polys::AbstractArray{T}) = (seriestype --> :shapefile; shapefile_coords(polys))
        @recipe f{T<:Shapefile.Handle}(handle::T) = handle.shapes
    end
end

# ----------------------------------------------------------------

# a recipe for composite lines with varying attributes
@userplot CompositeLine

@recipe function f(h::CompositeLine)
    x, y, linegroups = h.args

    seriesx = Dict{eltype(linegroups), Vector{Float64}}()
    seriesy = Dict{eltype(linegroups), Vector{Float64}}()

    for i in 1:length(linegroups)
        if ! haskey(seriesx, linegroups[i])
            seriesx[linegroups[i]] = Float64[]
            seriesy[linegroups[i]] = Float64[]
        end

        push!(seriesx[linegroups[i]], x[i])
        push!(seriesy[linegroups[i]], y[i])

        if i > 1 && linegroups[i] != linegroups[i-1]
            # End the previous series with the same point to avoid gaps
            push!(seriesx[linegroups[i-1]], x[i])
            push!(seriesy[linegroups[i-1]], y[i])

            # Add an NaN to stop the group
            push!(seriesx[linegroups[i-1]], NaN)
            push!(seriesy[linegroups[i-1]], NaN)

        end
    end

    seriestype := :path

    for key in keys(seriesx)
        @series begin
            seriesx[key], seriesy[key]
        end
    end

end


# -------------------------------------------------------------------
# AST trees

function add_ast(adjlist, names, depthdict, depthlists, nodetypes, ex::Expr, parent_idx)
    idx = length(names)+1
    iscall = ex.head == :call
    push!(names, iscall ? string(ex.args[1]) : string(ex.head))
    push!(nodetypes, iscall ? :call : :expr)
    l = Int[]
    push!(adjlist, l)

    depth = parent_idx==0 ? 1 : depthdict[parent_idx] + 1
    depthdict[idx] = depth
    while length(depthlists) < depth
        push!(depthlists, Int[])
    end
    push!(depthlists[depth], idx)

    for arg in (iscall ? ex.args[2:end] : ex.args)
        if isa(arg, Expr) && arg.head == :line
            continue
        end
        push!(l, add_ast(adjlist, names, depthdict, depthlists, nodetypes, arg, idx))
    end
    idx
end

function add_ast(adjlist, names, depthdict, depthlists, nodetypes, x, parent_idx)
    push!(names, string(x))
    push!(nodetypes, :leaf)
    push!(adjlist, Int[])
    idx = length(names)

    depth = parent_idx==0 ? 1 : depthdict[parent_idx] + 1
    depthdict[idx] = depth
    while length(depthlists) < depth
        push!(depthlists, Int[])
    end
    push!(depthlists[depth], idx)

    idx
end

@recipe function f(ex::Expr)
    names = String[]
    adjlist = Vector{Int}[]
    depthdict = Dict{Int,Int}()
    depthlists = Vector{Int}[]
    nodetypes = Symbol[]
    add_ast(adjlist, names, depthdict, depthlists, nodetypes, ex, 0)
    names := names
    # method := :tree
    method := :buchheim
    root --> :top

    markercolor --> Symbol[(nt == :call ? :pink : nt == :leaf ? :white : :lightgreen) for nt in nodetypes]

    # # compute the y-values from the depthdict dict
    # n = length(depthlists)-1
    # layers = Float64[(depthdict[i]-1)/n for i=1:length(names)]
    # # add_noise --> false
    #
    # positions = zeros(length(names))
    # for (depth, lst) in enumerate(depthlists)
    #     n = length(lst)
    #     pos = n > 1 ? linspace(0, 1, n) : [0.5]
    #     for (i, idx) in enumerate(lst)
    #         positions[idx] = pos[i]
    #     end
    # end
    #
    # layout_kw := KW(:layers => layers, :add_noise => false, :positions => positions)

    GraphPlot(get_source_destiny_weight(adjlist))
end


# -------------------------------------------------------------------
# Type trees

function add_subs!{T}(nodes, source, destiny, ::Type{T}, supidx)
    for sub in subtypes(T)
        push!(nodes, sub)
        subidx = length(nodes)
        push!(source, supidx)
        push!(destiny, subidx)
        add_subs!(nodes, source, destiny, sub, subidx)
    end
end

# recursively build a graph of subtypes of T
@recipe function f{T}(::Type{T}; namefunc = _->_.name.name)
    # get the supertypes
    sups = [T]
    sup = T
    while sup != Any
        sup = supertype(sup)
        unshift!(sups,sup)
    end

    # add the subtypes
    n = length(sups)
    nodes = copy(sups)
    source, destiny = collect(1:n-1), collect(2:n)
    add_subs!(nodes, source, destiny, T, n)

    # set up the graphplot
    names := map(namefunc, nodes)
    method --> :buchheim
    root --> :top
    GraphPlot((source, destiny))
end
