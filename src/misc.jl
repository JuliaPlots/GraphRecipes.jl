
# Remove legend, axes, grid, and force the aspect ratio for shape plots
@recipe function f(::Type{Val{:shapeplot}}, plt::Plot)
    legend --> false
    ticks --> nothing
    grid --> false
    aspect_ratio --> 1
    seriestype := :shape
    ()
end
@shorthands shapeplot

# ----------------------------------------------------------------
# Shapefile.jl plotting

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

        @recipe f(poly::Shapefile.ESRIShape) = (seriestype --> :shape; shapefile_coords(poly))
        @recipe f{T<:Shapefile.ESRIShape}(polys::AbstractArray{T}) = (seriestype --> :shape; shapefile_coords(polys))
    end
end

# ----------------------------------------------------------------
