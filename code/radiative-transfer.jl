include("common.jl")
m = KerrMetric(M = 1.0, a = 1.0)
x = SVector(0.0, 1000.0, deg2rad(80), 0.0)
d = PolishDoughnut(m)

Gradus.emissivity_coefficient(::AbstractMetric, ::PolishDoughnut, x, ν) = 0.1

pf = PointFunction((m, gp, t) -> gp.aux[1])

a, b, img = @time rendergeodesics(
    m,
    x,
    d,
    2000.0,
    verbose = true,
    pf = pf,
    αlims = (-25, 25), 
    βlims = (-15, 16),
    image_height = 800,
    image_width = 800,
    trace = Gradus.TraceRadiativeTransfer(I₀ = 0),
)

begin
    fig = Figure(size = (400, 300))
    ax = Axis(fig[1,1], aspect = DataAspect(), xlabel = "α", ylabel = "β", title = "Radiative transfer example")
    heatmap!(ax, a, b, img', colormap = :batlow)
    Makie.save("presentation/figs/raw/radiative-transfer.png", fig)
    fig
end