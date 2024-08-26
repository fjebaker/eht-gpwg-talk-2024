using Gradus, Plots

# specify the metric (AbstractMetric)
m = KerrMetric(a = 0.998)
m = JohannsenMetric(a = 0.6, ϵ3 = 2.0)
# Boyer-Lindquist coordinates for the observer
x_obs = SVector(0.0, 10_000.0, deg2rad(75), 0.0)
# choose a disc model (AbstractAccretionDisc)
d = ThinDisc(0.0, 20.0)

# compose a function that pulls out the quantities
# we are interested in visualizing
pf = ConstPointFunctions.redshift(m, x_obs) ∘
    # filter only those geodesics that intersected with
    # geometry
    ConstPointFunctions.filter_intersected()

# pass to a utility function, along with the maximum
# integration time
α, β, img = rendergeodesics(
  m, x_obs, d, 2x_obs[2]; 
  pf = pf, verbose = true,
  αlims = (-25, 25),
  βlims = (-20, 20),
)

# and visualise
Plots.heatmap(α, β, img, aspect_ratio=1, xlabel = "α", ylabel = "β")
Plots.savefig("presentation/figs/raw/code-example.svg")