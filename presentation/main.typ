#import "@preview/polylux:0.3.1": *
#import "tamburlaine.typ": *

#let HANDOUT_MODE = true
#enable-handout-mode(HANDOUT_MODE)

#show: tamburlaine-theme.with(aspect-ratio: "4-3")
#show link: item => underline(text(blue)[#item])

#let COLOR_CD = color.rgb("#56B4E9")
#let COLOR_REFL = color.rgb("#D55E00")
#let COLOR_CONT = color.rgb("#0072B2")

#set par(spacing: 0.5em, leading: 0.5em)
#set text(tracking: -0.5pt, size: 22pt)

#let uob_logo = read("./figs/UoB_CMYK_24.svg")

#let todo(item) = text(fill: red, [#item])

#title-slide(
  title_size: 30pt,
  title: [
    #set text(tracking: -3pt)
    #set align(left)
    #move(
      dy: -0.3em,
      text(
        size: 120pt,
        move(dy: -31pt, stack(
          spacing: 4pt,
          move(dx: -43pt, dy: 17pt, [#text(fill: TEXT_COLOR)[Spacetime]]),
          move(dx: 0pt, text(weight: "regular")[ agnostic]),
          [ray#text(weight: "regular")[_tracing_]]
        )
      ))
    )
    #v(0.5em)
  ],
  authors: ([Fergus Baker#super("1")], text(weight: "regular", [Andrew Young#super("1")])),
  where: "Grav. Phys. WG",
)[
  #align(right)[#image.decode(uob_logo, width: 20%)]
]

#slide(title: "Outline")[
  #set align(horizon)
  1. Gradus.jl
    - Spacetime agnostic approaches to general relativistic ray-tracing
    - Decomposing and piecing together model components
  2. Our immediate applications in X-ray astronomy
    - Relativistic blurring of reflection spectra
    - X-ray tests of GR?
  3. Possible applications for the GPWG
    - Catalogue of spacetimes
    - Expedite parameter exploration

  #v(1em)
]

#slide(title: "Gradus.jl")[
  #set align(horizon)
  #set align(center)
  #image(width: 30%, "./figs/gradus_logo.png")
  #v(1em)
  An open-source general relativistic ray-tracer written in Julia

  #v(1em)
  #align(right)[
    #set text(size: 15pt)
    Baker & Young, in prep.
  ]
]

#slide(title: "The geodesic equation")[
  // ray tracing all about solving geodesic equation
  #set align(horizon)

  #{
    set text(size: 30pt)
    $ dot.double(x)^mu + Gamma^(mu)_(nu sigma) dot(x)^nu dot(x)^sigma = 0 $
  }

  #v(1em)

  Initial value problem: set of 4 second order ODEs
  - Free choice of initial *three position* and *three velocity*
  - $dot(x)^t$ constrained by $g_(nu sigma) dot(x)^nu dot(x)^sigma = - mu^2$.

  #v(1em)

  Need *Christoffel components*
  - Want to use a *computer algebra system* to do this automatically
  - Simple for "simple" spacetimes
  - Can become horrendously complex for "non-simple" spacetimes
]

#slide(title: "Christoffel symbols on the fly")[
  #set align(horizon)

  $Gamma^mu_(nu sigma)$ is a function of $g_(mu nu)$ and the Jacobian w.r.t. the coordinates $partial_sigma [g_(mu nu)]$:
  - Use *automatic differentiation* to calculate $partial_sigma [g_(mu nu)]$ at every step
  - Fast and sufficiently accurate $cal(O)( epsilon N_"flop" )$

  #align(center)[
    #image(width: 50%, "./figs/christoffel-benchmark.svg")
  ]

  #v(1em)
  Flexibility: need only $g_(mu nu)$ rest we calculate numerically
  - General approach: implement all results derived from $g_(mu nu)$
]

#slide(title: "Toy accretion models")[
  #set align(horizon)

  First approximation of an *accretion disc* is the equatorial plane
  - Interested in set of (stable) time-like orbits in the equatorial plane
  - Have *analytic expressions* for $E$ and $L_z$ for certain spacetime classes (e.g. Johannsen, 2013)
  - Others can be solved with a mix of numerical methods, e.g. spacetimes immersed in a vector potential $A_mu$ (e.g. Hackmann et al., 2013)

  #v(1em)

  Can alternatively approach entirely numerically: orbit finding routines
  - Solve an optimization problem that minimizes $delta x^r$ and $delta x^theta$.
]

#slide(title: "Component modelling")[
  #set align(horizon)

  Design decision: maintain *public interfaces* for model components
  // what does this mean? it means you can pick and choose from a selection of
  // model components for building up your simulation
  // it means if you have a model for, e.g., the disc physics, you can just
  // focus on implementing that independent of the rest of the system
  // flexibility when exploring parameter spaces: one line changes to change the
  // entire simulation
  // separation of concerns: you can work on the analysis routines whilst
  // someone else works on a metric implementation, and the interfaces will
  // guarantee they work together
  // This can be done in any language, but is *easy* to do in Julia thanks to *multiple dispatch*

  #set text(size: 15pt)
  #grid(
    columns: (55%, 1fr),
    [
    ```julia
    using Gradus, Plots

    # specify the metric (AbstractMetric)
    m = KerrMetric(a = 0.998)
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
    heatmap(α, β, img, aspect_ratio=1)
    ```
    ],
    [
      #image(width: 80%, "./figs/code-example.svg")

      One line changes to update simulation:
      ```diff
      - m = KerrMetric(a = 0.998)
      + m = JohannsenMetric(a = 0.6, ϵ3 = 2.0)
      ```

      #image(width: 80%, "./figs/code-example-2.svg")
    ]
  )
  #v(1em)

]

#subtitle-slide[
  #set align(left)
  #text(size: 100pt)[
    #move(dy: -31pt, stack(
      spacing: 4pt,
      move(dx: -43pt, dy: 37pt, [#text(fill: TEXT_COLOR)[Our]]),
      move(dx: 0pt, text(weight: "regular")[_immediate_]),
      [#text(weight: "black")[applications]]
    ))
  ]
  #v(1em)
]

#slide(title: "Accretion in X-ray")[
  #v(-2.5em)
  #grid(
    columns: (58%, 40%),
    [
      #move(
        dx: -1em,
        animsvg(
          read("./figs/corona-plot.svg"),
          (i, im) => only(i)[
            #image.decode(im, width: 100%)
          ],
          (hide: ("g109", "g71", "g110", "g111")),
          (display: ("g110", "g111")),
          (hide: ("g110",)),
          (display: ("g109", "g71")),
          handout: HANDOUT_MODE,
        )
      )
    ],
    [
      #v(3.5em)
      Hello World
    ]
  )

  #set align(horizon)

  Measuring metric parameters using properties of the Fe K$alpha$ line
  - Reflection models given some illuminating flux

  Interested in how the *geometry* of the *(spacetime|disc|corona)* changes spectral profile

  Building models of *reflection spectra* and associated *reverberation lags*
]

#slide(title: "Our immediate needs")[
  #set align(horizon)

  Static, axisymmetric metrics of the form:
  $
  dif s^2 = g_(t t) dif t^2 + g_(r r) dif r^2 + g_(theta theta) dif theta^2 + g_(phi phi) dif phi^2  + 2 g_(t phi) dif t dif phi
  $

  We calculate *Cunningham transfer functions* to encode GR effects

  Involves solving for _specific_ geodesics, not just the image plane
  - Calculating the change-of-variable Jacobian

  Using *optimizers* to solve boundary value problems
  - *automatic differentiation* to trace dual numbers through the ODE system
]

#slide(title: "Applications on the horizon")[
  #set align(horizon)
  Things we are also looking to do
  - X-ray tests of relativity
  - Johannsen Psaltis metric
]

#subtitle-slide[
  #set align(left)
  #text(size: 100pt)[
    #move(dy: -31pt, stack(
      spacing: 4pt,
      move(dx: -43pt, dy: 37pt, [#text(fill: TEXT_COLOR)[Applications]]),
      move(dx: 0pt, text(weight: "regular")[_beyond_]),
      [#text(weight: "regular")[the] #text(weight: "black")[horizon]]
    ))
  ]
  #v(1em)
]

#slide(title: "Briefly looked at before")[
  #set align(horizon)

  Charged spacetimes

  Breaking some of our assumptions:
  - $bb(Z)_2$ metric

  Discontinuous spacetimes, Kerr with diffractive corona

  Only used coordinates that have singularities at the horizon
  - Use coordinate transforms that give horizon penetrating coordinates
  - Thinking about ways of making these transforms automatic
]

#slide(title: "GPWG")[
  #set align(horizon)

  Studying horizonless compact objects
  - Pattern speed (conroy 2023)
  - Photon rings / lensing bands
]

#slide(title: "Thank you")[
  #set align(horizon)
  Links, acknowledgements, etc.
]
