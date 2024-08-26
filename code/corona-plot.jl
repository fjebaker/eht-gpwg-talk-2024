import Random
include("common.jl")

function trace_corona_trajectories(
    m::AbstractMetric,
    d::AbstractAccretionGeometry,
    model;
    callback = domain_upper_hemisphere(),
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
    n_samples = 54,
    kwargs...,
)
    xs, vs, _ = Gradus.sample_position_direction_velocity(m, model, sampler, n_samples)
    tracegeodesics(
        m,
        xs,
        vs,
        d,
        1000;
        callback = callback,
        kwargs...,
    )
end

m = KerrMetric(1.0, 0.5)
model = LampPostModel(h = 10.0)
d = ShakuraSunyaev(m; eddington_ratio = 0.1)
dim = 26

alltraj = trace_corona_trajectories(m, d, model)

function is_intersected(sol) 
    intersected = (sol.prob.p.status[] == StatusCodes.IntersectedWithGeometry)
    if intersected
        return sol.u[end][2] < 21
    end
    false
end
traj = filter(is_intersected, alltraj.u) ; length(traj)
out_traj = filter(!is_intersected, alltraj.u) ; length(traj)

# lets get the endpoints, then find trajectories that trace from the observer to those points
endpoints = [
    SVector{4}(t.u[end][1:4]...) for t in traj
]

function trace_to_observer(target, m, x_obs, d)
    a, b, _, acc = Gradus.optimize_for_target(target, m, x_obs; callback = domain_upper_hemisphere(), β₀ = 0)
    @show acc
    v = map_impact_parameters(m, x_obs, a, b)
    tracegeodesics(m, x_obs, v, d, 200)
end

x_obs = SVector(0.0, 32, deg2rad(55), deg2rad(0))

sols = [trace_to_observer(t, m, x_obs, d) for t in endpoints] ;

begin
    fig = Figure(size = (800, 600))
    ax = Axis3(
        fig[1, 1],
        aspect = (1, 1, 1),
        limits = (-dim, dim, -dim, dim, -dim, dim),
        elevation = π / 23, #π / 12,
        azimuth = -deg2rad(65),
        viewmode = :fitzoom,
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        xlabelvisible = false,
        ylabelvisible = false,
        xspinewidth = 0,
        yspinewidth = 0,
        zspinewidth = 0,
    )
    hidedecorations!(ax)

    R = Gradus.inner_radius(m)
    bounding_sphere!(ax; R = R, color = :black)

    for r in range(Gradus.isco(m), 24.0, step = 3)
        plotring(ax, r; height = cross_section(d, r),  horizon_r = R, color = :black, dim = dim)
    end

    _palette = _default_palette()
    for sol in traj 
        plot_sol(
            ax,
            sol;
            color = popfirst!(_palette),
            horizon_r = R,
            dim = dim,
            show_intersect = false,
        ) 
    end
    for sol in out_traj 
        plot_sol(
            ax,
            sol;
            color = popfirst!(_palette),
            horizon_r = R,
            dim = dim,
            show_intersect = false,
        ) 
    end
    _palette = _default_palette()
    for sol in sols 
        plot_sol(
            ax,
            sol;
            color = popfirst!(_palette),
            horizon_r = R,
            dim = dim,
            show_intersect = false,
        ) 
    end

    Makie.save("presentation/figs/raw/corona-plot.pdf", fig)
    fig
end
