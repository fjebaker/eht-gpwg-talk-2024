using Revise
include("common.jl")
using Printf

function calc_lineprofile(m::AbstractMetric; inc = 60)
    x = SVector(0.0, 1000.0, deg2rad(inc), 0.0)
    d = ThinDisc(0.0, Inf)

    # maximal integration radius
    maxrₑ = 50.0

    # emissivity function
    ε(r) = r^(-3)

    # g grid to do flux integration over
    gs = range(0.0, 1.4, 500)
    _, flux = lineprofile(gs, ε, m, x, d, maxrₑ = maxrₑ, verbose = true)
    
    gs, flux
end

data = map([KerrMetric(a=0.7), JohannsenMetric(a=0.7, α13=1.0), JohannsenMetric(a=0.7, α13=-1.0)]) do m
    m, calc_lineprofile(m)
end

data2 = map([KerrMetric(a=0.7), JohannsenMetric(a=0.7, α13=1.0), JohannsenMetric(a=0.7, α13=-1.0)]) do m
    m, calc_lineprofile(m; inc = 30)
end

begin
    fig = Figure(size=(400, 350))
    ax = Axis(fig[1,1], xlabel = "E / E₀", ylabel = "Flux (arb.)", title= "Line profile for different spins")
    
    _palette = _default_palette()
    for info in data
        a, d = info
        c = popfirst!(_palette)
        lines!(ax, d..., color = c, linewidth = 2.5)
    end

    _palette = _default_palette()
    for info in data2
        a, d = info
        c = popfirst!(_palette)
        lines!(ax, d..., color = c, linewidth = 1.0)
    end
    
    xlims!(ax, 0, 1.4)
    # if i make this zero it seems to not draw certain parts of the figure 
    ylims!(ax, -1e-8, nothing)
    
    Makie.save("presentation/figs/raw/deformed-lineprofiles.svg", fig)
    fig
end

function a13_cond(M, a)
    -((M + sqrt(M^2 - a^2)) / M)^3
end

function calc_exclusion(as, ϵs)
    regions = zeros(Float64, (length(as), length(ϵs)))
    Threads.@threads for i in eachindex(as)
        a = as[i]
        for (j, ϵ) in enumerate(ϵs)
            m = JohannsenMetric(M = 1.0, a = a, α13 = ϵ)
            regions[i, j] = if a13_cond(m.M, m.a) > m.α13
                NaN
            else
                try
					Gradus.isco(m)
				catch
					-1
				end
            end
        end
    end
    regions
end

as = range(-1.0, 1.0, 200)
ϵs = range(-10, 10, 200)

img = @time calc_exclusion(as, ϵs)

curve = a13_cond.(1, as)

begin
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1,1], xlabel = "Spin a", ylabel = "Deformation Parameter α₁₃", title="ISCO and region of naked singularity")
    contourf!(
        ax,
        as, 
        ϵs, 
        img,
		colormap = :batlow,
		levels = collect(-1:13)
    )
    # lines!(as, curve, color = :black, linewidth = 3.0)
	xlims!(-1, 1)
	ylims!(-10, 10)
    Makie.save("presentation/figs/raw/deformed-isco.svg", fig)
    fig
end