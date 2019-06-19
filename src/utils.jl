
"""
This function builds a BezierCurve which leaves point p vertically upwards and
arrives point q vertically upwards. It may create a loop if necessary.
It assumes the view is [0,1]. That can be modified using the `xview` and
`yview` keyword arguments (default: `0:1`).
"""
function directed_curve(x1, x2, y1, y2; xview = 0:1, yview = 0:1, root::Symbol = :bottom)
    if root in (:left, :right)
        # flip x/y to simplify
        x1, x2, y1, y2, xview, yview = y1, y2, x1, x2, yview, xview
    end
    x = Float64[x1, x1]
    y = Float64[y1]

    minx, maxx = extrema(xview)
    miny, maxy = extrema(yview)
    dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    flip = root in (:top, :right)
    need_loop = (flip && y1 <= y2) || (!flip && y1 >= y2)

    # these points give the initial/final "rise"
    # note: this is a function of distance between points and axis scale
    y_offset = if need_loop
        0.3dist
    else
        min(0.3dist, 0.5*abs(y2-y1))
    end
    y_offset = max(0.02*(maxy-miny), y_offset)

    if flip
        # got the other direction
        y_offset *= -1
    end
    push!(y, y1 + y_offset)

    # try to figure out when to loop around vs just connecting straight
    if need_loop
        if abs(x2-x1) > 0.1 * (maxx - minx)
            # go between
            sgn = x2 > x1 ? 1 : -1
            x_offset = 0.5 * abs(x2-x1)
            append!(x, [x1 + sgn * x_offset, x2 - sgn * x_offset])
        else
            # add curve points which will create a loop
            x_offset = 0.3 * (maxx - minx) * (rand(Bool) ? 1 : -1)
            append!(x, [x1 + x_offset, x2 + x_offset])
        end
        append!(y, [y1 + y_offset, y2 - y_offset])
    end

    append!(x, [x2, x2])
    append!(y, [y2 - y_offset, y2])
    if root in (:left,:right)
        # flip x/y to simplify
        x,y = y,x
    end
    x, y
end


function shorten_segment(x1, y1, x2, y2, shorten)
    xshort = shorten * (x2-x1)
    yshort = shorten * (y2-y1)
    x1+xshort, y1+yshort, x2-xshort, y2-yshort
end


"""
Randomly pick a point to be the center control point of a bezier curve,
which is both equidistant between the endpoints and normally distributed
around the midpoint.
"""
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

function control_point(xi, xj, yi, yj, dist_from_mid)
    xmid = 0.5 * (xi+xj)
    ymid = 0.5 * (yi+yj)

    # get the angle of y relative to x
    theta = atan((yj-yi) / (xj-xi)) + 0.5pi

    # dist = sqrt((xj-xi)^2 + (yj-yi)^2)
    # dist_from_mid = curvature_scalar * 0.5dist

    # now we have polar coords, we can compute the position, adding to the midpoint
    (xmid + dist_from_mid * cos(theta),
     ymid + dist_from_mid * sin(theta))
end

# Function from Plots/src/components.jl
"get an array of tuples of points on a circle with radius `r`"
function partialcircle(start_θ, end_θ, n = 20, r = 1)
    Tuple{Float64,Float64}[(r*cos(u), r*sin(u)) for u in
                           range(start_θ, stop = end_θ, length = n)]
end

# for chord diagrams:
function arcshape(θ1, θ2)
    vcat(
        partialcircle(θ1, θ2, 15, 1.05),
        reverse(partialcircle(θ1, θ2, 15, 0.95))
    )
end

# x and y limits for arc diagram ()
function arcdiagram_limits(x, source, destiny)
    @assert length(x) >= 2
    margin = abs(0.1*(x[2] - x[1]))
    xmin, xmax = extrema(x)
    r = abs(0.5 * (xmax - xmin))
    mean_upside = mean(source .< destiny)
    ylims = if mean_upside == 1.0
        (-margin, r+margin)
    elseif mean_upside == 0.0
        (-r-margin, margin)
    else
        (-r-margin, r+margin)
    end
    (xmin-margin, xmax+margin), ylims
end

function islabel(item)
    ismissing(item) && return false
    (typeof(item) <: AbstractFloat && isnan(item)) && return false
    !in(item, (nothing, false, ""))
end
# From Plots/src/utils.jl
isnothing(x::Nothing) = true
isnothing(x) = false

# From Plots/src/Plots.jl
ignorenan_extrema(x) = Base.extrema(x)
# From Plots/src/utils.jl
ignorenan_extrema(x::AbstractArray{F}) where {F<:AbstractFloat} = NaNMath.extrema(x)
# From Plots/src/components.jl
function extrema_plus_buffer(v, buffmult = 0.2)
    vmin,vmax = extrema(v)
    vdiff = vmax-vmin
    buffer = vdiff * buffmult
    vmin - buffer, vmax + buffer
end
