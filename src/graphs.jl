
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

# ---------------------------------------------------------------------------
# Arc Diagram

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
