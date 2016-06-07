@userplot type MarginalHist
    args
end

@recipe function f(h::MarginalHist)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractVector) || !(typeof(h.args[2]) <: AbstractVector)
        error("Marginal Histograms should be given two vectors.  Got: $(typeof(h.args))")
    end
    
    # set up the subplots
    x, y = h.args
    legend := false
    link := :both
    ticks := [nothing :auto nothing]
    grid := false
    foreground_color_subplot := [RGBA(0,0,0,0) :match RGBA(0,0,0,0)]
    layout := @layout [tophist           _
                 	   hist2d{0.9w,0.9h} righthist]
    
    # main histogram2d
    @series begin
        seriestype := :histogram2d
        subplot := 2
        x, y
    end
    
    # these are common to both marginal histograms... borrow the first color from the fill gradient
    fillcolor := getColor(colorscheme(get(d, :fillcolor, Plots.default_gradient())))
    fillalpha := 0.3
    linealpha := 0.3
    
    # upper histogram
    @series begin
        seriestype := :histogram
        subplot := 1
        x
    end
    
    # right histogram
    @series begin
        seriestype := :histogram
        orientation := :h
        subplot := 3
        y
    end
end

# # now you can plot like:
# marginalhist(rand(1000), rand(1000))
