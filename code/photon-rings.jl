using Gradus
using Makie, CairoMakie

function ringfo(m)
    x = SVector(0.0, 1000.0, deg2rad(30), 0.0)
    pf = PointFunction( (m, gp, t) -> gp.aux.winding) ∘ 
        FilterPointFunction((m, gp, t) -> gp.aux.winding > 1, NaN)
    a,b,img = @time rendergeodesics(
        m,
        x,
        2000.0,
        αlims = (-8, 8),
        βlims = (-8, 8),
        trace = TraceWindings(),
        image_width = 800,
        image_height = 800,
        pf = pf
    )
end

m1 = JohannsenMetric(a = 0.9, α13 = 0.0)
m2 = JohannsenMetric(a = 0.9, α13 = 1.0)
m3 = JohannsenMetric(a = 0.9, α13 = -1.0)

d1 = ringfo(m1)
d2 = ringfo(m2)
d3 = ringfo(m3)

begin
    fig = Figure(size=(700,300))
    ax1 = Axis(fig[1,1], aspect = DataAspect(), xlabel = "α", ylabel = "β", title = "α₁₃ = 0.0")
    ax2 = Axis(fig[1,2], aspect = DataAspect(), xlabel = "α", title = "α₁₃ = 1.0")
    ax3 = Axis(fig[1,3], aspect = DataAspect(), xlabel = "α", title = "α₁₃ = -1.0")

    hideydecorations!(ax2, grid=false)
    hideydecorations!(ax3, grid=false)

    contourf!(ax1, d2[1], d2[2], d2[3]', colormap = :reds)
    contourf!(ax1, d3[1], d3[2], d3[3]', colormap = :greens)
    contourf!(ax1, d1[1], d1[2], d1[3]', colormap = :blues)
    dd = replace(d1[3]', NaN => 0)
    contour!(ax1, d1[1], d1[2], dd, color = :black)

    contourf!(ax2, d3[1], d3[2], d3[3]', colormap = :greens)
    contourf!(ax2, d1[1], d1[2], d1[3]', colormap = :blues)
    contourf!(ax2, d2[1], d2[2], d2[3]', colormap = :reds)
    dd = replace(d2[3]', NaN => 0)
    contour!(ax2, d2[1], d2[2], dd, color = :black)

    contourf!(ax3, d2[1], d2[2], d2[3]', colormap = :reds)
    contourf!(ax3, d1[1], d1[2], d1[3]', colormap = :blues)
    contourf!(ax3, d3[1], d3[2], d3[3]', colormap = :greens)
    dd = replace(d3[3]', NaN => 0)
    contour!(ax3, d3[1], d3[2], dd, color = :black)


    xlims!(ax1, -8.5, 8.5)
    ylims!(ax1, -8.5, 8.5)
    xlims!(ax2, -8.5, 8.5)
    ylims!(ax2, -8.5, 8.5)
    xlims!(ax3, -8.5, 8.5)
    ylims!(ax3, -8.5, 8.5)
    
    
    Makie.save("presentation/figs/raw/photon-rings.svg", fig)
    fig
end