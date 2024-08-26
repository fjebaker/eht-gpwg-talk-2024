using Gradus
using Makie, CairoMakie
using BenchmarkTools
using Statistics

function calc_times(m::AbstractMetric, xs, vs)
    # only the metric
    mtime = @benchmark Gradus.metric_components($m, xx) setup=(xx=SVector{2,Float64}(rand(xs)[2:3]...))
    # metric plus the jacobian simultaneously
    jtime =  @benchmark Gradus.metric_jacobian($m, xx) setup=(xx=SVector{2,Float64}(rand(xs)[2:3]...))
    # full geodesic equation
    gtime = @benchmark Gradus.geodesic_equation($m, xx, vv) setup=(xx=rand(xs); vv=rand(vs))
    (;mtime, jtime, gtime)  
end

m = KerrMetric(a = 0.998)

# generate some random positions and velocities for use in the benchmark only
# the precise values aren't important, so don't need to re-calcualted for each metric
xs = map(1:1000) do _
    SVector(rand(), 1.05 * Gradus.inner_radius(m) + rand(), π * rand(), 2π * rand())
end
vs = map(xs) do _
    SVector(rand(), rand(), rand(), rand())
end

kerr_times = calc_times(m, xs, vs)

m = JohannsenMetric(a = 0.6)
j1_times = calc_times(m, xs, vs)

m = JohannsenMetric(a = 0.6, α13 = 0.4)
j2_times = calc_times(m, xs, vs)

m = JohannsenMetric(a = 0.6, α13 = 0.4, ϵ3 = 1.0)
j3_times = calc_times(m, xs, vs)

m = BumblebeeMetric(a = 0.3, l = 2.0)
bee_times = calc_times(m, xs, vs)

data = (;
    x = [
        1, 1.8, 2.6, 3.3, 4.0
    ],
    y = [
        mean(kerr_times.gtime.times),
        mean(bee_times.gtime.times),
        mean(j1_times.gtime.times),
        mean(j2_times.gtime.times),
        mean(j3_times.gtime.times),
    ],
    y_metric = [
        mean(kerr_times.mtime.times),
        mean(bee_times.mtime.times),
        mean(j1_times.mtime.times),
        mean(j2_times.mtime.times),
        mean(j3_times.mtime.times),
    ],
    y_jac = [
        mean(kerr_times.jtime.times),
        mean(bee_times.jtime.times),
        mean(j1_times.jtime.times),
        mean(j2_times.jtime.times),
        mean(j3_times.jtime.times),
    ],
)

barkwargs = (;
    color = RGBAf(0.0, 0.0, 0.0, 0.0),
    strokecolor = :black,
    width = 0.2, 
    strokewidth = 1.0,
)

begin 
    fig = Figure(size=(500, 300))
    ax = Axis(fig[1,1],
        title = "Metric calculation times",
        ylabel = "CPU time (ns)",
        xticks = (data.x, ["Kerr", "Bumblebee", "Johannsen 0", "Johannsen 1", "Johannsen 2"]),
        xticklabelrotation = deg2rad(20)
    )
    ylims!(ax, 0, nothing)
    b1 = barplot!(
        ax,
        data.x .+ 0.2, data.y;
        barkwargs...
    )
    b2 = barplot!(
        ax,
        data.x, data.y_jac;
        barkwargs...
    )
    b3 = barplot!(
        ax,
        data.x .- 0.2, data.y_metric;
        barkwargs...
    )
    
    Legend(fig[1, 2], [b1, b2, b3], ["Metric", "Jac+Metric", "Geod. Eq."])

    fig
    Makie.save("presentation/figs/raw/christoffel-benchmark.svg", fig)
end
