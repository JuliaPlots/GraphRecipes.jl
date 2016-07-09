
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
