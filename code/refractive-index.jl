using Gradus, Makie, CairoMakie


m = KerrRefractive(a = 0.9, n = 0.95, corona_radius = 10.0)
x = SVector(0.0, 1000.0, deg2rad(70), 0.0)
d = ThinDisc(0.0, 18.0)

pf = PointFunction((m, gp, t) -> gp.x[2]) ∘ ConstPointFunctions.filter_intersected()
a, b, img = rendergeodesics(m, x, d, 2x[2], verbose = true, pf = pf, αlims = (-25, 25), βlims = (-15, 15), image_width = 1000, image_height = 1000)

m2 = KerrRefractive(a = 0.9, n = 1.05, corona_radius = 10.0)
a2, b2, img2 = rendergeodesics(m2, x, d, 2x[2], verbose = true, pf = pf, αlims = (-25, 25), βlims = (-15, 15), image_width = 1000, image_height = 1000)

begin
    fig = Figure(size = (350, 420), backgroundcolor = RGBAf(0.0,0.0,0.0,0.0))
    ga = fig[1,1] = GridLayout()
    ax = Axis(ga[1,1], aspect = DataAspect(), title = "Refractive index 0.95", xlabel = "α", ylabel = "β", backgroundcolor = RGBAf(0.0,0.0,0.0,0.0))
    heatmap!(ax, a, b, img', colormap = :batlow)
    hidexdecorations!(ax, grid = false)
    ylims!(ax, -15, 16)

    ax2 = Axis(ga[2,1], aspect = DataAspect(), title = "Refractive index 1.05", xlabel = "α", ylabel = "β", backgroundcolor = RGBAf(0.0,0.0,0.0,0.0))
    heatmap!(ax2, a, b, img2', colormap = :batlow)
    ylims!(ax2, -15, 16)
    
    rowgap!(ga, 5)
    Makie.save("presentation/figs/raw/refractive-index.png", fig, dpi = 400)
    fig
end