include("common.jl")
using Printf

function calc_lineprofile(m::AbstractMetric)
    x = SVector(0.0, 1000.0, deg2rad(50), 0.0)
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

data = map([0.0, 0.5, 0.75, 0.9, 0.998]) do a
    a, calc_lineprofile(KerrMetric(a=a))
end

begin
    fig = Figure(size=(400, 350))
    ax = Axis(fig[1,1], xlabel = "E / E₀", ylabel = "Flux (arb.)", title= "Line profile for different spins")
    
    _palette = _default_palette()
    for info in data
        a, d = info
        c = popfirst!(_palette)
        lines!(ax, d..., color = c)
        _text = Printf.@sprintf "a=%.1f" a
        text!(ax, 1.17, maximum(d[2]) * 0.95, text=_text, color = c)
    end
    
    xlims!(ax, 0, 1.4)
    # if i make this zero it seems to not draw certain parts of the figure 
    ylims!(ax, -1e-8, nothing)
    
    Makie.save("presentation/figs/raw/line-profile-kerr.svg", fig)
    fig
end