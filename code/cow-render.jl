using Gradus
using Makie, CairoMakie
using FileIO

using CoordinateTransformations, Rotations

teapot = load("data/utahteapot.stl")
mesh(teapot, color = :blue)

static_teapot = map(teapot) do triangle
    Tuple(0.5 .* SVector(p[1], p[2], p[3]) for p in triangle)
end

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1000.0, deg2rad(75), 0.0)

xfm = LinearMap(RotZ(0)) ∘ Translation(0.5, -6.0, -2.5) ∘ LinearMap(RotZ(-π/2))
apply_xfm(xfm, mesh) = map(tri -> xfm.(tri), mesh)
teapot_xfm = apply_xfm(xfm, static_teapot)
bbox = Gradus.bounding_box(teapot_xfm)
teapot_disc = MeshAccretionGeometry(teapot_xfm, bbox...)

a, b, cache = @time prerendergeodesics(
    m,
    x,
    teapot_disc,
    2000.0
    ;
    image_width = 300 * 4,
    image_height = 250 * 4,
    αlims = (-5, 13),
    βlims = (-7, 7),
    verbose = true,
)

pf = PointFunction((m, gp, t) -> gp.x[2]) ∘ ConstPointFunctions.filter_intersected()
img = apply(pf, cache)

pf2 = PointFunction((m, gp, t) -> gp.x[1]) ∘ FilterPointFunction((m, gp, t ) -> gp.status == StatusCodes.WithinInnerBoundary, NaN)
horizon_img = apply(pf2, cache)

begin
    fig = Figure(size = (400, 300), backgroundcolor = RGBAf(0.0,0.0,0.0,0.0))    
    ax = Axis(fig[1,1], xlabel = "α", ylabel = "β", aspect = DataAspect(), backgroundcolor = RGBAf(0.0,0.0,0.0,0.0))
    heatmap!(ax, a, b, horizon_img'; colormap = Reverse(:reds))
    heatmap!(ax, a, b, img'; colormap = :batlow)
    
    resize_to_layout!(fig)
    Makie.save("presentation/figs/raw/teapot.png", fig, dpi=400)
    fig
end