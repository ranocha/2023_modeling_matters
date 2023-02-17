### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c590d66e-8e24-4ecd-b6ef-67ec750624f2
begin
	# https://github.com/fonsp/Pluto.jl/wiki/%F0%9F%8E%81-Package-management#pattern-the-pkg-cell
	import Pkg
	redirect_stdio(stderr = devnull) do
    	Pkg.activate(@__DIR__)
	end
    Pkg.instantiate()
	using LinearAlgebra
	using OrdinaryDiffEq
	using ForwardDiff
	using PlutoUI
	using Plots
	plotlyjs() # use the plotlyjs backend of Plots.jl
	import PlotlyJS
end

# ╔═╡ 0029218c-50f3-44e3-b2e7-ef1e153d5798
md"""
# Modeling still matters: a surprising instance of catastrophic floating point errors in mathematical biology and numerical methods for ODEs

This Pluto notebook contains the numerical experiments reported in the paper
and all code required to reproduce numerical results and figures.

To reproduce the numerical experiments presented in this article, you need 
to install [Julia](https://julialang.org/). The numerical experiments presented 
in this article were performed using Julia v1.8.3.

Then, you need to start Julia in this directory and execute the following commands
in the REPL

```julia
import Pkg; Pkg.activate("."); Pkg.instantiate(); import Pluto; Pluto.run()
```

You can combine this with starting Julia from the command line as follows:

```bash
julia -e 'import Pkg; Pkg.activate("."); Pkg.instantiate(); import Pluto; Pluto.run()'
```

Then, the web server of [Pluto.jl](https://github.com/fonsp/Pluto.jl) should start
and open a browser window for you. There, you need to select the file `notebook.jl`.
Then, Pluto.jl should load the file and you should start to see some text. The
setup will take up to several minutes since Julia needs to install all dependencies
and execute the code. When everything is finished, the notebook is ready for interactive
use and exploration. In particular, you will be able to change several parameters
interactively and Julia will make sure that everything is sufficiently fast to be
updated on the fly.
"""

# ╔═╡ 11bc02fd-a89b-4187-b43e-9eb12fb2fc5b
md"""
## Modeling the inheritance of genes

This is the numerical solution of the original 3-component model

$$\begin{aligned}
    q_1' &= q_1^2 + q_1 q_2 + \frac{1}{4} q_2^2 - q_1, \\
    q_2' &= \frac{1}{2} q_2^2 + q_1 q_2 + 2 q_1 q_3 + q_2 q_3 - q_2, \\
    q_3' &= \frac{1}{4} q_2^2 + q_2 q_3 + q_3^2 - q_3,
\end{aligned}$$

using the fifth-order explicit Runge-Kutta method of Tsitouras
with absolute and relative tolerances $10^{-8}$.

If you run this notebook interactively, you can zoom in, move the plot,
hoover over the curves to get tooltips with the values you are pointing at
and so on.
"""

# ╔═╡ 087ceabc-9217-4355-8399-6d6630012890
md"""
The genotype proportions $q_i$ seem to converge to a steady state at first but
suddenly start to go to zero quite fast between $t = 30$ and
$t = 40$ --- although the sum of all components should be conserved 
(equal to one) for all times.
"""

# ╔═╡ ac2fb1b9-9aae-4838-be15-8a2cfb7ef5fb
md"""
Save plot above $(@bind save_fig_system3_original_Tsit5 CheckBox(default=false))
"""

# ╔═╡ 5aced9ae-749a-49cb-ac72-cb3971cd3981
md"""
The method of Tsitouras used above is the default explicit Runge-Kutta method
of [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and used
quite often in practice. Thus, one can expect that it works reliably. However,
maybe we have just found an instance where it fails. Thus, let's see whether
the gold-standard fifth-order explicit Runge-Kutta method of Dormand and Prince
works better.
"""

# ╔═╡ 714d0986-7fcf-40ed-838a-3fa3514ce2cb
md"""
The behavior is the same for short times. However, the numerical solution 
blows-up instead of going to zero between $t = 30$ and $t = 40$!
"""

# ╔═╡ e84abcc5-ef8d-4aad-89a3-a2041eab5535
md"""
Save plot above $(@bind save_fig_system3_original_DP5_1 CheckBox(default=false))
"""

# ╔═╡ f0619666-3e65-418c-81a5-6faeb86ce0c3
md"""
We can see this clearly if we zoom out a bit.
"""

# ╔═╡ e1df018b-e89e-4a94-8d2c-7e0460805f4f
md"""
Save plot above $(@bind save_fig_system3_original_DP5_2 CheckBox(default=false))
"""

# ╔═╡ cf8b4dbe-32eb-4fb1-8810-b1a1bbfc5958
md"""
## Dynamical System Analysis
"""

# ╔═╡ 95b844a8-5a70-4c1d-81fd-e984d12252ea
md"""
### Steady states of the original 3-component system

Here you can see a visualization of the non-negative steady states 
of the system. The plotting routine verifies that only the trivial steady
state $q = 0$ is stable.
"""

# ╔═╡ 04b73a56-e80a-419e-b4ff-abd9bdebf262
md"""
Save plot above $(@bind save_fig_system3_original_steadystates CheckBox(default=false))
"""

# ╔═╡ 9fdc84a7-3736-4214-8dc0-8201f365746b
md"""
### Steady states of the modified 3-component system

Here, we consider the steady states of the modified 3-component system

$$\begin{aligned}
    q_1' &= \frac{1}{4} q_2^2 - q_1 q_3, \\
    q_2' &= -\frac{1}{2} q_2^2 + 2 q_1 q_3, \\
    q_3' &= \frac{1}{4} q_2^2 - q_1 q_3,
\end{aligned}$$

derived in the paper. If you run this notebook interactively, you can
rotate the point of view and interact with the plot.
"""

# ╔═╡ b346bbb4-5bdc-4f6d-a538-af2b4442689c
md"""
Save plot above $(@bind save_fig_system3_modified_steadystates CheckBox(default=false))
"""

# ╔═╡ 3951294a-2fae-476e-a959-44a50761020c
md"""
### A reduced 2-component model

The reduced 2-component model

$$\begin{aligned}
    q_1' &= a q_1^2 + q_1 q_2 + (1-a) q_2^2 - q_1, \\
    q_2' &= (1-a) q_1^2 + q_1 q_2 + a q_2^2 - q_2,
\end{aligned}$$

with parameter $a = 0.7$ shows the same qualitative behavior 
as the original 3-component model in numerical experiments.
Again, we use absolute and relative tolerances $10^{-8}$.
"""

# ╔═╡ f47c7e86-e4c8-4641-80ba-57c8307c51d0
md"""
First, we apply the fifth-order method of Tsitouras.
"""

# ╔═╡ ad6af3b5-219e-4054-bc7a-4589a56da433
md"""
Save plot above $(@bind save_fig_system2_original_Tsit5 CheckBox(default=false))
"""

# ╔═╡ 141ecaa2-4c6e-4716-adc0-701281a0ee47
md"""
In this case, the method of Dormand and Prince also yields numerical solutions
going to zero.
"""

# ╔═╡ ee593071-d19a-4c9f-8b1e-dfdd1871dd41
md"""
However, there are many other methods with numerical solutions blowing-up again,
e.g., the sixth-order accurate method of Verner.
"""

# ╔═╡ 83eab3fe-74b1-4ced-8822-b641ec939432
md"""
Save plot above $(@bind save_fig_system2_original_Vern6 CheckBox(default=false))
"""

# ╔═╡ c4a8ea57-dd63-436e-a6ab-179000c7aca7
md"""
### Changing the floating point accuracy

To verify that floating point errors are indeed the cause of the catastrophic
long-term behavior of the original model, we change the floating point accuracy
from 64 bit (`Float64` in Julia) to other precisions.
"""

# ╔═╡ a4e4a6cf-12b9-4495-8f5f-feee5e117c51
md"""
#### `Float32`

We reduce the accuracy by using only 32 bit (`Float32` in Julia). 
Since the 32 bit machine precision is approximately 
$(round(eps(Float32), sigdigits=2)),
we increase the tolerance to $$10^{-7}$$.
"""

# ╔═╡ a355b03c-8513-40dd-a680-f1a23eafe810
md"""
First, we use the method of Tsitouras.
"""

# ╔═╡ 7ebcb1df-9cf7-48bb-8810-5d0d08aa3615
md"""
Indeed, the qualitative behavior is the same as before --- but the 
numerical solution leaves the steady state earlier than with 64 bit computations.
"""

# ╔═╡ 0e665b70-61c0-4245-ae8f-a2c1fe2034e3
md"""
Save plot above $(@bind save_fig_system3_original_Tsit5_Float32 CheckBox(default=false))
"""

# ╔═╡ 928f5e73-ecf3-4a1a-b17c-11074f8310d1
md"""
Next, we use the method of Dormand and Prince.
"""

# ╔═╡ 41793e67-1ba0-42df-8c2f-65e22b420fbb
md"""
Again, the qualitative behavior is the same as with `Float64` --- the
catastrophic blow-up just happens earlier.
"""

# ╔═╡ 4809cf73-2e1a-4785-bb5d-fd091d38a0c2
md"""
Save plot above $(@bind save_fig_system3_original_DP5_Float32 CheckBox(default=false))
"""

# ╔═╡ b5196d8f-04fb-437c-881e-0c5c5e7faa35
md"""
#### `BigFloat`

Now, we increase the accuracy to higher precision floating point numbers
`BigFloat` in Julia. For these experiments, we also use stricter tolerances
$$10^{-14}$$.
"""

# ╔═╡ dc1dd442-ec97-400e-88e9-c01b56156398
md"""
The method of Tsitouras does not change its qualitative behavior ---
the effect of floating point errors just becomes visible much later.
"""

# ╔═╡ 02927c54-d0ae-4bc8-afe0-7dc48753af8e
md"""
Save plot above $(@bind save_fig_system3_original_Tsit5_BigFloat CheckBox(default=false))
"""

# ╔═╡ 3fac63b1-11f3-4eff-a707-301a4a681e60
md"""
The same happens for the method of Dormand and Prince.
"""

# ╔═╡ e00eced4-465d-4c95-9744-fbaf53c41cf1
md"""
Save plot above $(@bind save_fig_system3_original_DP5_BigFloat CheckBox(default=false))
"""

# ╔═╡ 2d8f9ec7-af4e-4a7f-94c0-270614656090
md"""
#### Interactive exploration

You can also explore this phenomenon interactively.
You can change the tolerances `abstol`, `reltol`,
the floating point type, the algorithm and all other
options in the cell below and see the results.
"""

# ╔═╡ 3351f22b-e782-48cd-b12b-32a311b96aba
md"""
Here is the same setup for the 2-component model.
"""

# ╔═╡ b6076b6b-862a-419c-8c04-9e1f8c71f716
md"""
## Modified 2-component model

Here, we solve the modified 2-component model 

$$\begin{aligned}
    q_1' &= (1-a) \bigl( - q_1^2 + q_2^2\bigr), \\
    q_2' &= (1-a) \bigl( q_1^2 -  q_2^2 \bigr),
\end{aligned}$$

with the method of Tsitouras.
The numerical solutions remain stable, even with the larger default tolerances.
"""

# ╔═╡ 9ab0fbc8-3391-4832-850c-6ec536de34cd
md"""
Save plot above $(@bind save_fig_system2_modified_Tsit5 CheckBox(default=false))
"""

# ╔═╡ 6246fb20-f20c-4033-a202-8f348bdf1ea4
md"""
You can adapt the code below to explore other initial conditions and algorithms.
For this, you need to hoover over the the plot above or the space to its left.
Then, you should see an eye symbol --- to see the code, you need to click on
this symbol.
"""

# ╔═╡ 9db939bf-b6db-4cbd-9ead-72da978c2386
md"""
### Vector field of the original 2-component model

Here is a plot of the vector field of the right-hand side of the
original 2-component model.
"""

# ╔═╡ 8a04770b-57f2-4ad3-9015-89a196722474
md"""
We can see that the steady state $$(q_1^\star, q_2^\star) = (1/2, 1/2)$$
is attractive *as long as the solution stays in the manifold given by
$$q_1 + q_2 = 1$$*. Even a tiny perturbation of the total sum $$q_1 + q_2$$
moves the solution out of the invariant manifold and leads to an unstable behavior.
"""

# ╔═╡ c5c6adae-fc0f-4614-8676-ba7826846fae
md"""
Save plot above $(@bind save_fig_system2_original_vectorfield CheckBox(default=false))
"""

# ╔═╡ 0adab9e7-bea6-4379-8a0c-4c74632607d6
md"""
### Vector field of the modified 2-component model
"""

# ╔═╡ 139d3e12-449b-4e65-beda-565b24259de1
md"""
We can see that a perturbation of a steady state does not lead to unstable
behavior. Any perturbation perpendicular to the manifold of steady states
is in the direction of the invariant manifolds with fixed sum $q_1 + q_2$.
"""

# ╔═╡ b2643137-4a21-451a-bac9-7b0b213f7558
md"""
Save plot above $(@bind save_fig_system2_modified_vectorfield CheckBox(default=false))
"""

# ╔═╡ 55143b22-38e3-474b-bae3-67df17d9a6a2
md"""
## Modified 3-component model

The modified 3-component model 

$$\begin{aligned}
    q_1' &= \frac{1}{4} q_2^2 - q_1 q_3, \\
    q_2' &= -\frac{1}{2} q_2^2 + 2 q_1 q_3, \\
    q_3' &= \frac{1}{4} q_2^2 - q_1 q_3,
\end{aligned}$$

does also not result in any unexpected behavior anymore.
"""

# ╔═╡ d0ff9640-a9f4-4b0b-8011-550a9bc0d3a2
md"""
Save plot above $(@bind save_fig_system3_modified_Tsit5 CheckBox(default=false))
"""

# ╔═╡ 63fac493-7c86-4890-ac0c-a0424751a1ef
md"""
That's it - the end of the numerical experiments accompagnying our paper.
If you scroll further down below, you can find additional code used for the
visualizations shown above.
"""

# ╔═╡ 27dba1f7-2dbb-47e0-a1c0-540ab487cce8
md"""
# Appendix

If you are only interested in using the interactive visualizations in this
notebook, you do not need to read further. The following cells contain just
some source code used above.
"""

# ╔═╡ f8583865-1f19-430a-963f-4166cc405416
# v1.36.6 of Plots.jl can produce arrows with PlotlyJS but 
# v1.37 of plots.jl fails to do so, see
# https://github.com/JuliaPlots/Plots.jl/issues/4593
# Thus, we need to use an older version of Plots.jl and install
# it as above. However, the old version of Plots.jl yields errors
# in `Plots.plotlyjs_syncplot(plt)` so we need to copy their code
# and patch it with the hack `haskey(series_dict, :type) || continue`.
function __plotlyjs_syncplot(plt)
    plt[:overwrite_figure] && Plots.closeall()
    plt.o = PlotlyJS.plot()
    traces = PlotlyJS.GenericTrace[]
    for series_dict in Plots.plotly_series(plt)
		haskey(series_dict, :type) || continue
        plotly_type = pop!(series_dict, :type)
        series_dict[:transpose] = false
        push!(traces, PlotlyJS.GenericTrace(plotly_type; series_dict...))
    end
    PlotlyJS.addtraces!(plt.o, traces...)
    layout = Plots.plotly_layout(plt)
    w, h = plt[:size]
    PlotlyJS.relayout!(plt.o, layout, width = w, height = h)
    return plt.o
end

# ╔═╡ c7c9d0c9-1788-422d-991d-a228c0f84aa3
space = html"<br><br><br>";

# ╔═╡ 03cbabf2-10b9-4d2f-bc22-6c9b163b5f85
space

# ╔═╡ c526e16f-4f18-49d3-966f-e999f3e78dc4
space

# ╔═╡ 2faa807a-348a-4365-a104-75f72c433fde
space

# ╔═╡ 67fcf498-fcc3-46ea-a0f0-d54629a09fe3
# The original system with 3 components
# where the total sum is a second integral
function system3_original!(dq, q, parameters, t)
  q1, q2, q3 = q
  dq[1] = q1^2 + q1 * q2 + (1//4) * q2^2 - q1
  dq[2] = (1//2) * q2^2 + q1 * q2 + 2 * q3 * q1 + q2 * q3 - q2
  dq[3] = (1//4) * q2^2 + q2 * q3 + q3^2 - q3
  return nothing
end

# ╔═╡ 83ca4ed2-2d15-4b69-ab9e-ad2b5254c919
# For root-finding analysis etc.
system3_original!(dq, q) = system3_original!(dq, q, nothing, nothing)

# ╔═╡ b8e5ef20-2a48-457c-8119-093b52797f99
# The modified system with 3 components
# where the total sum is a first integral
function system3_modified!(dq, q, parameters, t)
  q1, q2, q3 = q
  dq[1] = (1//4) * q2^2 - q1 * q3
  dq[2] = (-1//2) * q2^2 + 2 * q1 * q3
  dq[3] = (1//4) * q2^2 - q1 * q3
  return nothing
end

# ╔═╡ 4ea31a9d-7ca7-44e1-9f3a-86c4542cb809
# For root-finding analysis etc.
system3_modified!(dq, q) = system3_modified!(dq, q, nothing, nothing)

# ╔═╡ c49ee69e-acbc-45d8-bb32-810f29ca6bc6
# The original system with 2 components
# where the total sum is a second integral
function system2_original!(dq, q, params, t)
  q1, q2 = q
  dq[1] = (7//10) * q1^2 + q1 * q2 + (3//10) * q2^2 - q1
  dq[2] = (7//10) * q2^2 + q1 * q2 + (3//10) * q1^2 - q2
  return nothing
end

# ╔═╡ 97c7e095-f667-4f25-be98-07fda1f8c7d1
# For root-finding analysis etc.
system2_original!(dq, q) = system2_original!(dq, q, nothing, nothing)

# ╔═╡ 6ccaf73e-d3ce-419e-b3e0-81de6c284359
# The modified system with 2 components
# where the total sum is a first integral
function system2_modified!(dq, q, params, t)
  q1, q2 = q
  dq[1] = (3//10) * (-q1^2 + q2^2)
  dq[2] = (3//10) * ( q1^2 - q2^2)
  return nothing
end

# ╔═╡ 8b221c92-ccf1-4066-b45a-2e76706e8f0b
# For root-finding analysis etc.
system2_modified!(dq, q) = system2_modified!(dq, q, nothing, nothing)

# ╔═╡ b66eacce-9cc7-42b6-a9f5-f84fdb8907b0
plot_kwargs = let
	# fontsizes = (
	# 	xtickfontsize = 18, ytickfontsize = 18, 
	# 	xguidefontsize = 20, yguidefontsize = 20, 
	# 	legendfontsize = 18)
	fontsizes = (
		xtickfontsize = 14, ytickfontsize = 14, 
		xguidefontsize = 16, yguidefontsize = 16, 
		legendfontsize = 14)
	(; linewidth = 3, gridlinewidth = 2, fontsizes...)
end

# ╔═╡ f9c0e157-cee9-418c-9c9e-65c9df6157cf
plot_3d_kwargs = let
	fontsizes = (
		xtickfontsize = 9, ytickfontsize = 9, ztickfontsize = 9,
		xguidefontsize = 11, yguidefontsize = 11, zguidefontsize = 11,
		legendfontsize = 11)
	(; linewidth = 5, gridlinewidth = 3, fontsizes...)
end

# ╔═╡ 91b29fe5-0550-4deb-b9da-9881e67df0c9
let
	# ^2 to get a better spacing of the points
	q1 = range(0.0, 1.0, length = 201).^2
	q2 = range(0.0, 1.0, length = 201).^2
	q3 = @. q2^2 / (4 * q1')
	
	fig = plot(xguide = "q₁", yguide = "q₂", zguide = "q₃")
	plot!(fig, q1, q2, q3; seriestype = :surface,
		  clims = (0.0, 1.0), colorbar_title = "2-parameter family",
		  plot_3d_kwargs...)
	plot!(fig, zero(q1), zero(q1), q1; label = "q₁ = q₂ = 0", plot_3d_kwargs...,
		  linewidth = 10, color = palette(:default)[1])
	plot!(zrange = (0.0, 1.0))

	# steady states with ∑ᵢ qᵢ = 1
	q1 = range(0.0, 1.0, length = 201)
	q2 = @. 2 * (sqrt(q1) - q1)
	q3 = @. 1 - q1 - q2
	plot!(fig, q1, q2, q3; label="q₁ + q₂ + q₃ = 1", plot_3d_kwargs...,
		  linewidth = 12)
	
	plot!(fig, camera = (45, 50))
	global fig_system3_modified_steadystates = fig
end

# ╔═╡ 7e61e991-6125-4a25-8550-541a4f80f297
let
	if save_fig_system3_modified_steadystates
		local filename = savefig(
			fig_system3_modified_steadystates, 
			"fig_system3_modified_steadystates.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ caec4b21-9448-4390-a9e1-f59f77dcb019
function plot_single_sol(sol)
	t = range(sol.prob.tspan..., length = 200)
	fig = plot(; xguide = "Time t", yguide = "Genotype proportions qᵢ")
	plot!(fig, t, sol.(t, idxs = 1); label = "q₁", linestyle = :solid, plot_kwargs...)
	plot!(fig, t, sol.(t, idxs = 2); label = "q₂", linestyle = :dash, plot_kwargs...)
	if length(sol.u[1]) > 2
		plot!(fig, t, sol.(t, idxs = 3); label = "q₃", linestyle = :dashdot, plot_kwargs...)
	end
	return fig
end

# ╔═╡ f31150a3-81cb-4c12-9670-b2c513b5de3c
let
	q0 = [0.5, 0.25, 0.25]
	tspan = (0, 50)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, Tsit5(), abstol = 1.0e-8, reltol = 1.0e-8)
	global fig_system3_original_Tsit5 = plot_single_sol(sol)
end

# ╔═╡ b0a3abfd-dcf1-4a22-add1-13ec666f4910
let
	if save_fig_system3_original_Tsit5
		filename = savefig(
			fig_system3_original_Tsit5, "fig_system3_original_Tsit5.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 4fffe198-6e1b-4f03-90eb-be61b026b26c
let
	q0 = [0.5, 0.25, 0.25]
	tspan = (0, 50)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, DP5(), abstol = 1.0e-8, reltol = 1.0e-8)
	fig = plot_single_sol(sol)
	plot!(fig, ylims = (0.0, 1.0))
	global fig_system3_original_DP5 = fig
end

# ╔═╡ f20b8c7d-4d40-4c7b-9138-05a4869960f6
let
	if save_fig_system3_original_DP5_1
		local filename = savefig(
			fig_system3_original_DP5, "fig_system3_original_DP5_1.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 2a89f261-be0f-43f4-8797-4e0c77a3ad6b
plot(fig_system3_original_DP5, yscale = :log10, ylims = :auto, legend = :right)

# ╔═╡ be3a2c41-cd92-41fd-87dd-26604b21795a
let
	if save_fig_system3_original_DP5_2
		local filename = savefig(
			plot(fig_system3_original_DP5, yscale = :log10, ylims = :auto, legend = :right), 
			"fig_system3_original_DP5_2.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 3572f11d-8631-473c-9a5d-b6fa46037b14
let
	q0 = [0.25, 0.75]
	tspan = (0, 50)
	ode = ODEProblem(system2_original!, q0, tspan)
	sol = solve(ode, Tsit5(), abstol = 1.0e-8, reltol = 1.0e-8)
	fig = plot_single_sol(sol)
	global fig_system2_original_Tsit5 = fig
end

# ╔═╡ a2e3f854-513d-456e-ad36-49a8467c92e6
let
	if save_fig_system2_original_Tsit5
		local filename = savefig(
			fig_system2_original_Tsit5, 
			"fig_system2_original_Tsit5.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ d57f5714-73f6-45e7-9684-fa053ac966c0
let
	q0 = [0.25, 0.75]
	tspan = (0, 50)
	ode = ODEProblem(system2_original!, q0, tspan)
	sol = solve(ode, DP5(), abstol = 1.0e-8, reltol = 1.0e-8)
	fig = plot_single_sol(sol)
end

# ╔═╡ 3a0b5670-71ba-45fd-a7b3-993dc65ae8aa
let
	q0 = [0.25, 0.75]
	tspan = (0, 50)
	ode = ODEProblem(system2_original!, q0, tspan)
	sol = solve(ode, Vern6(), abstol = 1.0e-8, reltol = 1.0e-8)
	fig = plot_single_sol(sol)
	plot!(fig, ylim = (0.0, 1.0))
	global fig_system2_original_Vern6 = fig
end

# ╔═╡ 2d6bddb5-79bf-4370-8ddb-6643c0070c22
let
	if save_fig_system2_original_Vern6
		local filename = savefig(
			fig_system2_original_Vern6, 
			"fig_system2_original_Vern6.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 38782a5f-a506-420a-8343-5b233770e17f
let
	q0 = Float32[0.5, 0.25, 0.25]
	tspan = (0, 50)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, Tsit5(), abstol = 1.0e-7, reltol = 1.0e-7)
	fig = plot_single_sol(sol)
	global fig_system3_original_Tsit5_Float32 = fig
end

# ╔═╡ f5aa3cdf-4cc0-4590-86b7-972b425b611f
let
	if save_fig_system3_original_Tsit5_Float32
		local filename = savefig(
			fig_system3_original_Tsit5_Float32, 
			"fig_system3_original_Tsit5_Float32.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ c2c44e8a-596a-4a69-8eb5-bf23ed2b1fce
let
	q0 = Float32[0.5, 0.25, 0.25]
	tspan = (0, 50)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, DP5(), abstol = 1.0e-7, reltol = 1.0e-7)
	fig = plot_single_sol(sol)
	plot!(fig, ylims = (0.0, 1.0))
	global fig_system3_original_DP5_Float32 = fig
end

# ╔═╡ 6395f6c7-f9e5-49fc-b47d-781d4add4b7c
let
	if save_fig_system3_original_DP5_Float32
		local filename = savefig(
			fig_system3_original_DP5_Float32, 
			"fig_system3_original_DP5_Float32.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 2ecc495f-b4ab-4984-ba37-4490870a5367
let
	q0 = BigFloat[0.5, 0.25, 0.25]
	tspan = (0, 250)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, Tsit5(), abstol = 1.0e-14, reltol = 1.0e-14)
	fig = plot_single_sol(sol)
	global fig_system3_original_Tsit5_BigFloat = fig
end

# ╔═╡ 9b941007-7d5f-4569-91ee-b8e8e873d3ba
let
	if save_fig_system3_original_Tsit5_BigFloat
		local filename = savefig(
			fig_system3_original_Tsit5_BigFloat, 
			"fig_system3_original_Tsit5_BigFloat.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 4d19d588-8a33-4270-9aae-138df904ac07
let
	q0 = BigFloat[0.5, 0.25, 0.25]
	tspan = (0, 250)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, DP5(), abstol = 1.0e-14, reltol = 1.0e-14)
	fig = plot_single_sol(sol)
	plot!(fig, ylims = (0.0, 1.0))
	global fig_system3_original_DP5_BigFloat = fig
end

# ╔═╡ 58a339c8-6ced-4db8-bba9-6b8107ba0400
let
	if save_fig_system3_original_DP5_BigFloat
		local filename = savefig(
			fig_system3_original_DP5_BigFloat, 
			"fig_system3_original_DP5_BigFloat.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 180030fd-037b-4576-aa98-ab5586762566
let
	q0 = Float64[0.5, 0.25, 0.25] # Float32, Float64, BigFloat
	tspan = (0, 250)              # or just (0, 50)
	ode = ODEProblem(system3_original!, q0, tspan)
	sol = solve(ode, 
		Rodas5(); # Tsit5(), DP5(), ...
		abstol = 1.0e-12, reltol = 1.0e-12)
	fig = plot_single_sol(sol)
	# plot!(fig, ylims = (0.0, 1.0)) # if the solutions blows up
end

# ╔═╡ c71e24e2-67e1-4541-8142-b0948ce99627
let
	q0 = Float64[0.25, 0.75]
	tspan = (0, 250)
	ode = ODEProblem(system2_original!, q0, tspan)
	sol = solve(ode, Tsit5(), abstol = 1.0e-8, reltol = 1.0e-8)
	fig = plot_single_sol(sol)
	# plot!(fig, ylims = (0.0, 1.0)) # if the solutions blows up
end

# ╔═╡ 23e40a4c-8606-4196-ad32-a38c5615ce43
let
	q0 = [0.25, 0.75] # initial condition
	tspan = (0, 50)   # time span of the integraiton
	ode = ODEProblem(system2_modified!, q0, tspan)
	sol = solve(ode, Tsit5())
	fig = plot_single_sol(sol)
	global fig_system2_modified_Tsit5 = fig
end

# ╔═╡ e2ddfaff-7549-40ad-a949-16291f145d94
let
	if save_fig_system2_modified_Tsit5
		local filename = savefig(
			fig_system2_modified_Tsit5, 
			"fig_system2_modified_Tsit5.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 543bd7f2-1d63-4b23-9c6f-ae215f904a67
let
	q0 = [0.5, 0.25, 0.25]
	tspan = (0, 50)
	ode = ODEProblem(system3_modified!, q0, tspan)
	sol = solve(ode, Tsit5())
	fig = plot_single_sol(sol)
	plot!(fig, legend = :right)
	global fig_system3_modified_Tsit5 = fig
end

# ╔═╡ c003f4d7-daa7-45d5-b000-5a51e6579ba6
let
	if save_fig_system3_modified_Tsit5
		local filename = savefig(
			fig_system3_modified_Tsit5, 
			"fig_system3_modified_Tsit5.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 5a822c5d-2c5f-4739-9031-350f7949f655
function stable(max_real_λ; tol = sqrt(eps(eltype(max_real_λ))))
  if max_real_λ > tol
    return "unstable"
  else
    return "stable"
  end
end

# ╔═╡ 129b34ba-040c-455d-84f5-5bf5006ccde0
let
	q1 = range(0.0, 1.0, length = 31)
	q1 = q1.^2 # to spread the points more evenly
	q2 = @. 2 * (sqrt(q1) - q1)
	q3 = @. 1 + q1 - 2 * sqrt(q1)

	q0 = [q1[1], q2[1], q3[1]]
	cfg = ForwardDiff.JacobianConfig(system3_original!, similar(q0), q0)
	J = similar(q0, length(q0), length(q0))
  	j! = (J, q) -> ForwardDiff.jacobian!(J, system3_original!, q0, q, cfg)
	max_real_λ = similar(q1)
	for (i, q) in enumerate(Iterators.zip(q1, q2, q3))
		@. q0 = q
		j!(J, q0)
		λ = eigvals(J)
		max_real_λ[i] = maximum(real, λ)
	end
	
	fig = plot(xguide = "q₁", yguide = "q₂", zguide = "q₃")
	plot!(fig, q1, q2, q3; label="unstable", plot_3d_kwargs...)
	if all(==("unstable"), stable.(max_real_λ))
		println("all states are unstable")
	end
	scatter!(fig, [0], [0], [0]; label="stable", plot_3d_kwargs...)
	plot!(fig, camera = (55, 10), legend = :right, bottom_margin = -10 * Plots.mm)
	global fig_system3_original_steadystates = fig
end

# ╔═╡ 13a5acd3-924b-4fa2-a2f7-7cffc40cf45c
let
	if save_fig_system3_original_steadystates
		local filename = savefig(
			fig_system3_original_steadystates, 
			"fig_system3_original_steadystates.pdf")
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ 0103c06a-a071-48a2-afd6-090530a6707d
meshgrid(x, y) = (repeat(x, outer = length(y)), repeat(y, inner = length(x)))

# ╔═╡ cd3f7ad3-1f92-40b1-9d8d-23bf4c6b284f
let
	scale!(f) = (f .*= 0.06 ./ sqrt(norm(f)))
	n = 21
	q1 = range(0.0, 1.0, length = n)
	q2 = range(0.0, 1.0, length = n)
	q1, q2 = meshgrid(q1, q2)
	f1 = zero(q1)
	f2 = zero(q2)
	f = zeros(2)
	for (i, q) in enumerate(Iterators.zip(q1, q2))
		system2_original!(f, q, nothing, nothing)
		scale!(f)
		f1[i], f2[i] = f
	end
	color = range(0.0, 1.0, length = length(f1))
	fig = plot(xguide = "q₁", yguide = "q₂", aspect_ratio = :equal,
			   size = (600, 600))
	plot!([1.0, 0.0], [0.0, 1.0]; label = "q₁ + q₂ = 1", color = 2,
		  plot_kwargs...)
	quiver!(fig, q1, q2, quiver = (f1, f2);
			plot_kwargs..., linewidth = 2, color = 1, linecolor = 1)
	# Linecolors are sadly not supported by the PlotlyJS backend...
	# fig = quiver(q1, q2, quiver = (f1, f2), line_z = color)
	global fig_system2_original_vectorfield = fig
end

# ╔═╡ 00f7c0c4-fee8-4f8e-b42d-80d796d1fd85
let
	if save_fig_system2_original_vectorfield
		# local filename = savefig(
		# 	fig_system2_original_vectorfield, 
		# 	"fig_system2_original_vectorfield.pdf")
		# Currently, there is a bug in Plots.jl such that we
		# need the following workaround, see
		# https://github.com/JuliaPlots/Plots.jl/issues/3199
		plt = fig_system2_original_vectorfield
		width, height = plt.attr[:size]
		Plots.prepare_output(plt)
		local filename = PlotlyJS.savefig(
			# This does not work with Plots.jl v1.36.6, 
			# but v1.37 does not work correctly with arrows, see
			# https://github.com/JuliaPlots/Plots.jl/issues/4593
			# Plots.plotlyjs_syncplot(plt), 
			__plotlyjs_syncplot(plt), 
			"fig_system2_original_vectorfield.pdf", 
			width=width, height=height)
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ e1862809-7960-4dd6-a146-1ae4c3ccbb66
let
	# gr()
	scale!(f) = (f .*= 0.06 ./ sqrt(norm(f)))
	n = 21
	q1 = range(0.0, 1.0, length = n)
	q2 = range(0.0, 1.0, length = n)
	q1, q2 = meshgrid(q1, q2)
	f1 = zero(q1)
	f2 = zero(q2)
	f = zeros(2)
	for (i, q) in enumerate(Iterators.zip(q1, q2))
		system2_modified!(f, q, nothing, nothing)
		scale!(f)
		f1[i], f2[i] = f
	end
	color = range(0.0, 1.0, length = length(f1))
	fig = plot(xguide = "q₁", yguide = "q₂", aspect_ratio = :equal,
			   size = (600, 600))
	plot!(fig, [1.0, 0.0], [0.0, 1.0]; label = "q₁ + q₂ = 1", color = 2,
		  plot_kwargs...)
	quiver!(fig, q1, q2, quiver = (f1, f2);
			plot_kwargs..., linewidth = 2, color = 1, linecolor = 1)
	plot!(fig, [0.0, 1.0], [0.0, 1.0]; label = "steady states", color = 3,
		  plot_kwargs...)
	# Linecolors are sadly not supported by the PlotlyJS backend...
	# fig = quiver(q1, q2, quiver = (f1, f2), line_z = color)
	# plotlyjs()
	global fig_system2_modified_vectorfield = fig
end

# ╔═╡ 2bd8020b-d028-4cbf-ae94-b3530707645e
let
	if save_fig_system2_modified_vectorfield
		# local filename = savefig(
		# 	fig_system2_modified_vectorfield, 
		# 	"fig_system2_modified_vectorfield.pdf")
		# Currently, there is a bug in Plots.jl such that we
		# need the following workaround, see
		# https://github.com/JuliaPlots/Plots.jl/issues/3199
		plt = fig_system2_modified_vectorfield
		width, height = plt.attr[:size]
		Plots.prepare_output(plt)
		local filename = PlotlyJS.savefig(
			# This does not work with Plots.jl v1.36.6, 
			# but v1.37 does not work correctly with arrows, see
			# https://github.com/JuliaPlots/Plots.jl/issues/4593
			# Plots.plotlyjs_syncplot(plt), 
			__plotlyjs_syncplot(plt), 
			"fig_system2_modified_vectorfield.pdf", 
			width=width, height=height)
		Markdown.parse("""Plot saved as `$(filename)` """)
	end
end

# ╔═╡ Cell order:
# ╟─0029218c-50f3-44e3-b2e7-ef1e153d5798
# ╟─11bc02fd-a89b-4187-b43e-9eb12fb2fc5b
# ╟─f31150a3-81cb-4c12-9670-b2c513b5de3c
# ╟─087ceabc-9217-4355-8399-6d6630012890
# ╟─ac2fb1b9-9aae-4838-be15-8a2cfb7ef5fb
# ╟─b0a3abfd-dcf1-4a22-add1-13ec666f4910
# ╟─5aced9ae-749a-49cb-ac72-cb3971cd3981
# ╟─4fffe198-6e1b-4f03-90eb-be61b026b26c
# ╟─714d0986-7fcf-40ed-838a-3fa3514ce2cb
# ╟─e84abcc5-ef8d-4aad-89a3-a2041eab5535
# ╟─f20b8c7d-4d40-4c7b-9138-05a4869960f6
# ╟─f0619666-3e65-418c-81a5-6faeb86ce0c3
# ╟─2a89f261-be0f-43f4-8797-4e0c77a3ad6b
# ╟─e1df018b-e89e-4a94-8d2c-7e0460805f4f
# ╟─be3a2c41-cd92-41fd-87dd-26604b21795a
# ╟─cf8b4dbe-32eb-4fb1-8810-b1a1bbfc5958
# ╟─95b844a8-5a70-4c1d-81fd-e984d12252ea
# ╟─129b34ba-040c-455d-84f5-5bf5006ccde0
# ╟─04b73a56-e80a-419e-b4ff-abd9bdebf262
# ╟─13a5acd3-924b-4fa2-a2f7-7cffc40cf45c
# ╟─9fdc84a7-3736-4214-8dc0-8201f365746b
# ╟─91b29fe5-0550-4deb-b9da-9881e67df0c9
# ╟─b346bbb4-5bdc-4f6d-a538-af2b4442689c
# ╟─7e61e991-6125-4a25-8550-541a4f80f297
# ╟─3951294a-2fae-476e-a959-44a50761020c
# ╟─f47c7e86-e4c8-4641-80ba-57c8307c51d0
# ╟─3572f11d-8631-473c-9a5d-b6fa46037b14
# ╟─ad6af3b5-219e-4054-bc7a-4589a56da433
# ╟─a2e3f854-513d-456e-ad36-49a8467c92e6
# ╟─141ecaa2-4c6e-4716-adc0-701281a0ee47
# ╟─d57f5714-73f6-45e7-9684-fa053ac966c0
# ╟─ee593071-d19a-4c9f-8b1e-dfdd1871dd41
# ╟─3a0b5670-71ba-45fd-a7b3-993dc65ae8aa
# ╟─83eab3fe-74b1-4ced-8822-b641ec939432
# ╟─2d6bddb5-79bf-4370-8ddb-6643c0070c22
# ╟─c4a8ea57-dd63-436e-a6ab-179000c7aca7
# ╟─a4e4a6cf-12b9-4495-8f5f-feee5e117c51
# ╟─a355b03c-8513-40dd-a680-f1a23eafe810
# ╟─38782a5f-a506-420a-8343-5b233770e17f
# ╟─7ebcb1df-9cf7-48bb-8810-5d0d08aa3615
# ╟─0e665b70-61c0-4245-ae8f-a2c1fe2034e3
# ╟─f5aa3cdf-4cc0-4590-86b7-972b425b611f
# ╟─928f5e73-ecf3-4a1a-b17c-11074f8310d1
# ╟─c2c44e8a-596a-4a69-8eb5-bf23ed2b1fce
# ╟─41793e67-1ba0-42df-8c2f-65e22b420fbb
# ╟─4809cf73-2e1a-4785-bb5d-fd091d38a0c2
# ╟─6395f6c7-f9e5-49fc-b47d-781d4add4b7c
# ╟─b5196d8f-04fb-437c-881e-0c5c5e7faa35
# ╟─dc1dd442-ec97-400e-88e9-c01b56156398
# ╟─2ecc495f-b4ab-4984-ba37-4490870a5367
# ╟─02927c54-d0ae-4bc8-afe0-7dc48753af8e
# ╟─9b941007-7d5f-4569-91ee-b8e8e873d3ba
# ╟─3fac63b1-11f3-4eff-a707-301a4a681e60
# ╟─4d19d588-8a33-4270-9aae-138df904ac07
# ╟─e00eced4-465d-4c95-9744-fbaf53c41cf1
# ╟─58a339c8-6ced-4db8-bba9-6b8107ba0400
# ╟─2d8f9ec7-af4e-4a7f-94c0-270614656090
# ╠═180030fd-037b-4576-aa98-ab5586762566
# ╟─3351f22b-e782-48cd-b12b-32a311b96aba
# ╠═c71e24e2-67e1-4541-8142-b0948ce99627
# ╟─b6076b6b-862a-419c-8c04-9e1f8c71f716
# ╟─23e40a4c-8606-4196-ad32-a38c5615ce43
# ╟─9ab0fbc8-3391-4832-850c-6ec536de34cd
# ╟─e2ddfaff-7549-40ad-a949-16291f145d94
# ╟─6246fb20-f20c-4033-a202-8f348bdf1ea4
# ╟─9db939bf-b6db-4cbd-9ead-72da978c2386
# ╟─cd3f7ad3-1f92-40b1-9d8d-23bf4c6b284f
# ╟─8a04770b-57f2-4ad3-9015-89a196722474
# ╟─c5c6adae-fc0f-4614-8676-ba7826846fae
# ╟─00f7c0c4-fee8-4f8e-b42d-80d796d1fd85
# ╟─0adab9e7-bea6-4379-8a0c-4c74632607d6
# ╟─e1862809-7960-4dd6-a146-1ae4c3ccbb66
# ╟─139d3e12-449b-4e65-beda-565b24259de1
# ╟─b2643137-4a21-451a-bac9-7b0b213f7558
# ╟─2bd8020b-d028-4cbf-ae94-b3530707645e
# ╟─55143b22-38e3-474b-bae3-67df17d9a6a2
# ╟─543bd7f2-1d63-4b23-9c6f-ae215f904a67
# ╟─d0ff9640-a9f4-4b0b-8011-550a9bc0d3a2
# ╟─c003f4d7-daa7-45d5-b000-5a51e6579ba6
# ╟─63fac493-7c86-4890-ac0c-a0424751a1ef
# ╟─03cbabf2-10b9-4d2f-bc22-6c9b163b5f85
# ╟─c526e16f-4f18-49d3-966f-e999f3e78dc4
# ╟─2faa807a-348a-4365-a104-75f72c433fde
# ╟─27dba1f7-2dbb-47e0-a1c0-540ab487cce8
# ╠═c590d66e-8e24-4ecd-b6ef-67ec750624f2
# ╠═f8583865-1f19-430a-963f-4166cc405416
# ╠═c7c9d0c9-1788-422d-991d-a228c0f84aa3
# ╠═67fcf498-fcc3-46ea-a0f0-d54629a09fe3
# ╠═83ca4ed2-2d15-4b69-ab9e-ad2b5254c919
# ╠═b8e5ef20-2a48-457c-8119-093b52797f99
# ╠═4ea31a9d-7ca7-44e1-9f3a-86c4542cb809
# ╠═c49ee69e-acbc-45d8-bb32-810f29ca6bc6
# ╠═97c7e095-f667-4f25-be98-07fda1f8c7d1
# ╠═6ccaf73e-d3ce-419e-b3e0-81de6c284359
# ╠═8b221c92-ccf1-4066-b45a-2e76706e8f0b
# ╠═b66eacce-9cc7-42b6-a9f5-f84fdb8907b0
# ╠═f9c0e157-cee9-418c-9c9e-65c9df6157cf
# ╠═caec4b21-9448-4390-a9e1-f59f77dcb019
# ╠═5a822c5d-2c5f-4739-9031-350f7949f655
# ╠═0103c06a-a071-48a2-afd6-090530a6707d
