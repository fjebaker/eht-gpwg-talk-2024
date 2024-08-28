using Gradus
using Printf, Makie, CairoMakie

function NoZDisc(m::AbstractMetric{T}; outer_radius = T(50)) where {T}
    isco = Gradus.isco(m)

    rs = collect(range(isco, outer_radius, 200))
    hs = rs .* cos.(Gradus._solve_orbit_θ.(m, rs))

    interp = Gradus._make_interpolation(rs, hs)
    WarpedThinDisc(interp, isco, outer_radius)
end

function apparent_image(m)
    x = SVector(0.0, 10000.0, deg2rad(85), 0.0)
    d = NoZDisc(m, outer_radius = 10.0)
    pf = ConstPointFunctions.redshift(m, x) ∘ ConstPointFunctions.filter_intersected()
    α, β, img = rendergeodesics(
        m,
        x,
        d,
        # maximum integration time
        20000.0,
        αlims = (-13, 13),
        βlims = (-10, 10),
        image_width = 800,
        image_height = 800,
        verbose = true,
        pf = pf,
    )
end

m1 = NoZMetric(ϵ = 0.0, a = 0.9)
m2 = NoZMetric(ϵ = 2.0, a = 0.9)
m3 = NoZMetric(ϵ = -2.0, a = 0.9)
m4 = NoZMetric(ϵ = 1.0, a = 0.9)

d1 = apparent_image(m1)
d2 = apparent_image(m2)
d3 = apparent_image(m3)
d4 = apparent_image(m4)


begin
    fig = Figure(size=(800, 300))
    ax1 = Axis(fig[1,1], aspect=DataAspect(), title = "a=0.9, ε=0.0", xlabel = "α", ylabel = "β")
    ax2 = Axis(fig[1,2], aspect=DataAspect(), title = "a=0.9, ε=2.0", xlabel = "α")
    ax3 = Axis(fig[1,3], aspect=DataAspect(), title = "a=0.9, ε=-2.0", xlabel = "α")
    ax4 = Axis(fig[1,4], aspect=DataAspect(), title = "a=0.9, ε=1.0", xlabel = "α")
    
    for (i, ax) in enumerate((ax1, ax2, ax3, ax4))
        xlims!(ax, -13, 13)  
        ylims!(ax, -9.5, 9.5)  
        if i != 1
            hideydecorations!(ax, grid=false)
        end
    end
    
    ckwargs = (;
        levels = 50, 
        colormap = :batlow
    )
    
    contourf!(ax1, d1[1], d1[2], d1[3]'; ckwargs...)
    contourf!(ax2, d2[1], d2[2], d2[3]'; ckwargs...)
    contourf!(ax3, d3[1], d3[2], d3[3]'; ckwargs...)
    contourf!(ax4, d4[1], d4[2], d4[3]'; ckwargs...)

    Makie.save("presentation/figs/raw/noz-metric-images.svg", fig)
    fig
end
