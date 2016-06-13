@userplot CorrPlot

@recipe function f(cp::CorrPlot)
	mat = cp.args[1]
	n = size(mat,2)
	C = cor(mat)

	seriestype := :scatter
	link := :both
	layout := (n,n)
	t = vcat(:auto, fill(nothing, n-1))
	xticks := vec(repmat(reverse(t), 1, n))'
	yticks := t'

	sts = fill(:scatter, n, n)
	xs = Array(Any, n, n)
	ys = Array(Any, n, n)
	for i=1:n
		vi = mat[:,i]

		# histograms on the diagonal
		@series begin
			seriestype := :histogram
			subplot := ((i-1)*n + 1)
			vi
		end

		for j = i+1:n
			vj = mat[:,j]

			# under the diagonal... scatter
			@series begin
				subplot := ((i-1)*n + j)
				vj, vi
			end
		end
	end
end