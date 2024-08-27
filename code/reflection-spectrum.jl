using Reflionx
using SpectralFitting

include("common.jl")

reftable = Reflionx.parse_run("data/reflionx/grid")

labels = [
    "logξ = 1.5",
    "logξ = 2.0",
    "logξ = 2.5"
]

begin
    fig = Figure(size=(500, 350))  
    ax = Axis(fig[1,1], yscale = log10, ylabel = "Flux (arb.)", xlabel = "Energy (keV)", title = "Reflection spectra for different emissivities")

    vspan!(ax, [6.2], [6.5], color = (:grey, 0.3))

    datas = map(1:3) do i
        ref_spec = reftable.grids[5, i, 1]
        # upscale the grid we are interested in
        interp = DataInterpolations.LinearInterpolation(ref_spec.flux, ref_spec.energy)

        erange = collect(range(0.01, 10.1, 901))
        values = interp.(erange)[1:end-1] ./ diff(erange)
        erange = erange[1:end-1]
        lines!(ax, erange, values, label = labels[i])
    end
    
    xlims!(ax, 0, 10)
    axislegend(ax)
    
    Makie.save("presentation/figs/raw/reflection-spectrum.svg", fig)
    fig
end