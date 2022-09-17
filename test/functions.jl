function random_labelled_graph()
    n = 15
    seed!(RNG, 1)
    A = Float64[ rand(RNG) < 0.5 ? 0 : rand(RNG) for i=1:n, j=1:n]
    for i=1:n
        A[i, 1:i-1] = A[1:i-1, i]
        A[i, i] = 0
    end
    x = rand(RNG,n)
    y = rand(RNG,n)
    z = rand(RNG,n)
    p = graphplot(A,
              nodesize = 0.2,
              node_weights = 1:n,
              nodecolor = range(colorant"yellow", stop=colorant"red", length=n),
              names = 1:n,
              fontsize = 10,
              linecolor = :darkgrey,
              layout_kw = Dict(:x => x, :y => y),
              rng=RNG
              )
    p, n, A, x, y, z
end

function random_3d_graph()
    seed!(RNG, 1)
    n, A, x, y, z = random_labelled_graph()[2:end]
    graphplot(A,
              node_weights = 1:n,
              markercolor = :darkgray,
              dim = 3,
              markersize = 5,
              markershape = :circle,
              linecolor = :darkgrey,
              linealpha = 0.5,
              layout_kw = Dict(:x => x, :y => y, :z => z),
              rng=RNG,
              )
end

function light_graphs()
    seed!(RNG, 1)
    g = wheel_graph(10)
    graphplot(g, curves=false, rng=RNG)
end

function directed()
    g = [0 1 1;
         0 0 1;
         0 1 0]

    graphplot(g, names=1:3, curvature_scalar=0.1, rng=RNG)
end

function edgelabel()
    seed!(RNG, 1)
    n = 8
    g = wheel_digraph(n)
    edgelabel_dict = Dict()
    for i in 1:n
        for j in 1:n
            edgelabel_dict[(i, j)] = string("edge ", i, " to ", j)
        end
    end

    graphplot(g, names=1:n, edgelabel=edgelabel_dict, curves=false, nodeshape=:rect, rng=RNG)
end

function selfedges()
    seed!(RNG, 1)
    g = [1 1 1;
         0 0 1;
         0 0 1]

    graphplot(DiGraph(g), self_edge_size=0.2, rng=RNG)
end

function multigraphs()
    seed!(RNG, 1)
    graphplot([[1,1,2,2],[1,1,1],[1]], names="node_".*string.(1:3), nodeshape=:circle, self_edge_size=0.25, rng=RNG)
end

function arc_chord_diagrams()
    seed!(RNG, 2)
    adjmat = Symmetric(sparse(rand(RNG,0:1,8,8)))
    plot(
         graphplot(adjmat,
                   method=:chorddiagram,
                   names=[text(string(i), 8) for i in 1:8],
                   linecolor=:black,
                   fillcolor=:lightgray,
                   rng=RNG),

         graphplot(adjmat,
                   method=:arcdiagram,
                   markersize=0.5,
                   markershape=:circle,
                   linecolor=:black,
                   markercolor=:black,
                   rng=RNG)
    )
end

function marker_properties()
    N = 8
    seed = 42
    seed!(RNG, seed)
    g = barabasi_albert(N, 1; rng=RNG, seed=seed)
    weights = [length(neighbors(g, i)) for i in 1:nv(g)]
    graphplot(g, curvature_scalar=0,
              node_weights=weights, nodesize=0.25,
              linecolor=:gray,
              linewidth=2.5,
              nodeshape=:circle,
              node_z=rand(RNG,N), markercolor=:viridis,
              nodestrokewidth=1.5,
              markerstrokestyle=:solid,
              markerstrokealpha=1.0,
              markerstrokecolor=:lightgray,
              colorbar=true,
              rng=RNG,
    )
end

function ast_example()
    code = :(
    function mysum(list)
        out = 0
        for value in list
            out += value
        end
        out
    end
    )

    plot(code, fontsize=10, shorten=0.01, axis_buffer=0.15, nodeshape=:rect, size=(1000, 1000), rng=RNG)
end

function julia_type_tree()
    plot(AbstractFloat, method=:tree, fontsize=10, nodeshape=:ellipse, size=(1000, 1000), rng=RNG)
end

AbstractTrees.children(d::Dict) = [p for p in d]
AbstractTrees.children(p::Pair) = AbstractTrees.children(p[2])
function AbstractTrees.printnode(io::IO, p::Pair)
    str = isempty(AbstractTrees.children(p[2])) ? string(p[1], ": ", p[2]) : string(p[1], ": ")
    print(io, str)
end

function julia_dict_tree()
    d = Dict(:a => 2,:d => Dict(:b => 4,:c => "Hello"),:e => 5.0)

    plot(TreePlot(d), method=:tree, fontsize=10, nodeshape=:ellipse, size=(1000, 1000), rng=RNG)
end

function diamond_nodeshape(x_i, y_i, s)
    [(x_i + 0.5s*dx, y_i + 0.5s*dy) for (dx, dy) in [(1,1),(-1,1),(-1,-1),(1,-1)]]
end
function diamond_nodeshape_wh(x_i, y_i, h, w)
    out = Tuple{Float64, Float64}[(-0.5,0),(0,-0.5),(0.5,0),(0,0.5)]
    map(out) do t
        x = t[1]* h 
        y = t[2]* w
        (
         x + x_i, 
         y + y_i 
        )
    end
end

function custom_nodeshapes_single()
    seed!(RNG, 6)
    g = rand(RNG,5,5)
    g[g .> 0.5] .= 0
    for i in 1:5
        g[i,i] = 0
    end
    graphplot(g, nodeshape=diamond_nodeshape, rng=RNG)
end

function custom_nodeshapes_various()
    seed!(RNG, 6)
    g = rand(RNG,5,5)
    g[g .> 0.5] .= 0
    for i in 1:5
        g[i,i] = 0
    end
    graphplot(g, nodeshape=[:circle, diamond_nodeshape, diamond_nodeshape_wh, :hexagon, diamond_nodeshape_wh], rng=RNG)
end

function funky_edge_and_marker_args()
    n = 5
    g = SimpleDiGraph(n)

    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 4)
    add_edge!(g, 4, 5)

    curviness_matrix = zeros(n,n)
    edgewidth_matrix = zeros(n,n)
    edgestyle_dict = Dict()
    for e in edges(g)
        curviness_matrix[e.src,e.dst] = 0.5sin(e.src)
        edgewidth_matrix[e.src,e.dst] = 0.8e.dst
        edgestyle_dict[(e.src,e.dst)] = e.src < 2.0 ? :solid :
                                        e.src > 3.0 ? :dash : :dot
    end
    edge_z_function = (s,d,w) -> curviness_matrix[s,d]

    graphplot(g, names=[" I "," am "," a ","funky","graph"],
        x=[1,2,3,4,5],
        y=[5,4,3,2,1],
        nodesize=0.3,
        size=(1000,1000),
        axis_buffer=0.6,
        fontsize=16,
        self_edge_size=1.3,
        curvature_scalar=curviness_matrix,
        edgestyle=edgestyle_dict,
        edgewidth=edgewidth_matrix,
        edge_z=edge_z_function,
        nodeshape=:circle,
        node_z=[1,2,3,4,5],
        nodestroke_z=[5,4,3,2,1],
        edgecolor=:viridis,
        markercolor=:viridis,
        nodestrokestyle=[:dash,:solid,:dot],
        nodestrokewidth=6,
        linewidth=2,
        colorbar=true,
        rng=RNG,
    )
end
