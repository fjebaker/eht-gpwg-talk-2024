#import "@preview/polylux:0.3.1": *
#import "tamburlaine.typ": *

#let HANDOUT_MODE = false
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
#let cbox(content, ..args) = rect(radius: 3pt, outset: 5pt, ..args, content)

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

// #slide(title: "Outline")[
//   #set align(horizon)
//   1. Gradus.jl
//     - Spacetime agnostic approaches to general relativistic ray-tracing
//     - Decomposing and piecing together model components
//   2. Our immediate applications
//     - Accretion in X-ray: relativistic reflection spectra
//     - Our needs and solutions
//     - X-ray tests of GR?
//   3. Possible applications for the GPWG
//     - Catalogue of spacetimes
//     - Expedite parameter exploration

//   #v(1em)
// ]

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

  #uncover("2-")[
    Need *Christoffel components*
    - Use a *computer algebra system* to do this automatically?
    #uncover("3-")[
    - Can be simple for "simple" spacetimes
    - Can become horrendous for "non-simple" spacetimes
    ]
  ]
]

#slide(title: "Christoffel symbols on the fly")[
  #v(0.5em)

  $Gamma^mu_(nu sigma)$ is a function of $g_(mu nu)$ and the Jacobian w.r.t. the coordinates $partial_sigma [g_(mu nu)]$:
  - Use *automatic differentiation* to calculate $partial_sigma [g_(mu nu)]$ at every step
  - Fast and sufficiently accurate $cal(O)( epsilon N_"flop" )$

  #grid(
    columns: (50%, 1fr),
    [
    #align(center)[
      #image(width: 100%, "./figs/christoffel-benchmark.svg")
    ]
    *Flexibility*: need only $g_(mu nu)$ rest we calculate numerically
    ],
    align(center)[
      #image(width: 90%, "./figs/deflection-angle.svg")
    ]
  )

  #v(0.5em)
  #cbox(fill: PRIMARY_COLOR, width: 90%, text(fill: SECONDARY_COLOR)[
    - Generalized approach: implement *all algorithms* directly from $g_(mu nu)$.
  ])

]

#slide(title: "Example: toy accretion models")[
  #v(1em)
  First approximation of an *accretion disc* is the equatorial plane
  #set text(size: 20pt)
  - Interested in set of (stable) *closed time-like orbits* in the equatorial plane
  - Have *analytic expressions* for $E$ and $L_z$ for certain spacetime classes (e.g. Johannsen, 2013)
  #v(0.5em)
  #grid(
    columns: (50%, 1fr),
    [
      #v(0.5em)
      Can alternatively approach *entirely numerically*: orbit finding routines
    - Solve an *optimization problem* that minimizes $delta x^r$ and $delta x^theta$.
    ]
  )

  #v(1fr)
  Others classes of spacetimes can be solved with a *mix of numerical methods* (Baker & Young, in prep.)
  - E.g. singularity immersed in a vector potential $A_mu$ (e.g. Hackmann et al., 2013)
  #v(2em)
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
    #one-by-one(start: 2, [
    ```julia
    using Gradus, Plots

    # specify the metric (AbstractMetric)
    m = KerrMetric(a = 0.998)
    # Boyer-Lindquist coordinates for the observer
    x_obs = SVector(0.0, 10_000.0, deg2rad(75), 0.0)
    # choose a disc model (AbstractAccretionDisc)
    d = ThinDisc(0.0, 20.0)
    ```
    #v(1em)
  ],[
    ```julia
    # compose a function that pulls out the quantities
    # we are interested in visualizing
    pf = ConstPointFunctions.redshift(m, x_obs) ∘
        # filter only those geodesics that intersected with
        # geometry
        ConstPointFunctions.filter_intersected()
    ```
    #v(1em)
  ],[
    ```julia
    # pass to a utility function, along with the maximum
    # integration time
    α, β, img = rendergeodesics(
      m, x_obs, d, 2x_obs[2];
      pf = pf, verbose = true,
      αlims = (-25, 25),
      βlims = (-20, 20),
    )
    ```
    #v(1em)
  ],[
    ```julia
    # and visualise
    heatmap(α, β, img, aspect_ratio=1)
    ```
  ])
    ],
    [
      #uncover("5-")[
      #image(width: 80%, "./figs/code-example.svg")
      ]

      #uncover("6-")[
      One line changes to update simulation:
      ```diff
      - m = KerrMetric(a = 0.998)
      + m = JohannsenMetric(a = 0.6, ϵ3 = 2.0)
      ```

      #image(width: 80%, "./figs/code-example-2.svg")
      ]
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
    columns: (58%, 1fr),
    column-gutter: 10pt,
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
          (),
          (),
          handout: HANDOUT_MODE,
        )
      )
      #v(1em)
      #set text(size: 20pt)
      #uncover("4-")[
      Convolve *line profile* with the reflection spectrum to embed GR effects
      ]
      #v(1em)
      #uncover("5-")[
      Measuring *metric parameters* using properties of the Fe K$alpha$ line (6.4 keV)
      #v(0.5em)
      ]
      #uncover("6-")[
      #cbox(fill: PRIMARY_COLOR, width: 95%, text(fill: SECONDARY_COLOR)[
      - Interested in how *geometry* of *(spacetime|disc|corona)* changes spectra
      - Relativistic models of *reflection spectra* and associated *reverberation lags*
      ])
      ]
    ],
    [
      #v(3.0em)
      #uncover("3-")[
      #image(width: 90%, "./figs/reflection-spectrum.svg")
      #text(size: 12pt)[Reflection spectra using version of `reflionx` of Ross, 1978.]
      ]

      #uncover("4-")[
      #image(width: 85%, "./figs/line-profile-kerr.svg")
      ]
    ]
  )
]

#slide(title: "Use case: transfer functions")[
  #set text(size: 18pt)
  #v(1em)
  #grid(
    columns: (55%, 1fr),
    column-gutter: 10pt,
    [
      #animsvg(
        read("./figs/ip-parameterization.svg"),
        (i, im) => only(i)[
          #image.decode(im, width: 90%)
        ],
        (hide: ("g184",)),
        (display: ("g184",)),
        (),
        (),
        (),
        (),
        (),
        handout: HANDOUT_MODE,
      )
      #uncover("2-")[
      Encode *relativistic effects* in *transfer functions* (Cunningham, 1975):
      re-parameterize image plane:
      $ (alpha, beta) arrow.r (r_"em", g^star) $
      ]

      #uncover("5-")[
      - Can be pre-computed, efficiently interpolated
      ]

      #uncover("6-")[
      #v(0.5em)
      Transfer functions can be *augmented* in various ways (Baker & Young, in prep):
      ]
      #uncover("7-")[
      - Include *timing* information, *optical depth* from radiative transfer, *polarization*, image parameters, and so on...
      ]
    ],
    [
      $
      F(E_"obs") &= integral.double I(E_"obs") dif alpha dif beta \
      #uncover("3-")[
      $
      arrow.b \
      F(E_"obs") &= integral.double I(E_"obs") #text(fill: PRIMARY_COLOR)[$abs( (partial (alpha, beta)) / (partial (r_"em", g^star)) )$] dif r_"em" dif g^star$ \
    ]
      $
      #v(1em)
      #uncover("4-")[
      With *Liouville's Theorem*, can now entirely express physics in the *frame of the disc*, and use transfer functions to "ray-trace back".
      ]

      #uncover("7-")[
      #image("./figs/radiative-transfer.svg", width: 90%)
      ]
    ]
  )
]

#slide(title: "What we needed")[
  #set align(horizon)

  *Static, axisymmetric* spacetimes with metrics of the form:
  $
  dif s^2 = g_(t t) dif t^2 + g_(r r) dif r^2 + g_(theta theta) dif theta^2 + g_(phi phi) dif phi^2  + 2 g_(t phi) dif t dif phi
  $

  #v(1em)
  #uncover("2-")[
  We mainly concern ourselves with *optically thick* accretion discs.
  ]


  #v(1em)

  #uncover("3-")[
  For transfer functions: interested in solving for *specific geodesics*, not full image plane
  - Using *optimizers* to solve boundary value problems (e.g. mapping impact parameters)
  - *Automatic differentiation* to trace dual numbers along the geodesics (e.g. Jacobian terms)
  ]

  #v(1em)

  #uncover("4-")[
  Good *verification schemes*: multiple approaches to the same problem
  - Implement *conceptually simple*, but *computational slow* methods to check the *conceptually complex*, but *computationally fast*
  - Extensive and ever-growing *test-suite*
  ]
]

#slide(title: "Applications on the horizon")[
  #set text(size: 20pt)

  #grid(
    columns: (50%, 1fr),
    [
      X-ray tests of GR:
      - Using *reflection spectra* to measure metric parameters
      #align(center)[
        // TODO: this needs a legend
        #image("./figs/deformed-lineprofiles.svg", width: 70%)
      ]

      #v(1em)
      #uncover("3-")[
      Line profiles can become *degenerate* when ISCOs are equal
      ]
      #uncover("4-")[
        #cbox(width: 90%, fill: PRIMARY_COLOR, text(fill: SECONDARY_COLOR)[
      - Can *timing* lift that degeneracy?
      - Do more complete models lift the degeneracy?
        ])
      ]
    ],
    [
      #set align(center)
      #uncover("2-")[
      #image("./figs/bambi-figure.svg", width: 80%)
      #align(right)[#text(size: 12pt)[Above figure from Bambi et al., 2023 #h(2em)]]
      ]
      #uncover("3-")[
        #image("./figs/deformed-isco.svg", width: 90%)
      ]
    ],
  )


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

#slide(title: "Examined before")[
  #set text(size: 20pt)
  As part of *developing* and *testing* Gradus.jl:
  - Explored charged spacetimes (e.g. Kerr-Newman, $Q eq.not 0$)
  - Photon rings and lensing bands within other spacetimes:
  #grid(
    column-gutter: 20pt,
    columns: (60%, 1fr),
    [
      #align(center)[
        #image("./figs/photon-rings.svg", width: 100%)
      ]
      #uncover("3-")[
      Only used *coordinates* that have *singularities along horizon*:
      - Not a concrete assumption in the implementation
      ]
      #uncover("4-")[
      - There are *coordinate transforms* that give horizon penetrating coordinates
      - Thinking about methods for making these transforms *automatic*
      ]
    ],
    [
      #set text(size: 18pt)
      #v(0.5em)
      #uncover("2-")[
      Discontinuous spacetimes, e.g. Kerr with refractive corona:

      #image("./figs/refractive-index.png", width: 95%)
      ]
    ]
  )


]

#slide(title: "Breaking symmetries")[
  No $bb(Z)_2$ metric:
  #align(center)[
    #image("./figs/chen-pu-24.png", width: 70%)
    #align(right, text(size: 12pt)[Above figure from Chen & Pu, 2024])
    #move(dx:-26pt, image("./figs/noz-metric-images.svg", width: 70%))
  ]
  #set text(size: 18pt)
  Useful for *validating other GRRT* codes
  #uncover("2-")[
  - Export transfer function tables and fit models to data directly #uncover("3-")[(_well_, almost!)]
  ]
  #align(center)[
    #uncover("4-")[
      #image("./figs/github-bug-report.png", width: 57%)
    ]
  ]
]

#slide(title: "Catalogue of spacetimes")[
  #uncover("6-")[
  Maintain a library of spacetimes:
  ]

  #grid(
    columns: (50%, 1fr),
    column-gutter: 10pt,
    row-gutter: 0.5em,
    [
      #uncover("6-")[
      - Minkowski
      - Kerr
      - Kerr + Dark Matter
      ]
    ],
    [
      #uncover("6-")[
      - Kerr + Refractive Bubble
      - Kerr-Newman
      - Morris-Thorne Wormhole
      ]
    ],
    [
      #v(0.5em)
      Implement a *new metric*:
      #v(0.5em)
      #set text(size: 15pt)
      #one-by-one(start: 2, [
      ```julia
      # `T` lets us pick number type
      struct Schwarzschild{T} <:
              AbstractStaticAxisSymmetric{T}
          M::T
      end
      ```
      #v(1em)
    ],
    [
      ```julia
      # event horizon; boundary of integration
      Gradus.inner_radius(m::Schwarzschild) = 2 * m.M
      ```
      #v(1em)
    ],
    [
      ```julia
      # define the static, axis-symmetric components
      function Gradus.metric_components(m::Schwarzschild, x)
          r, θ = x
          dt = -(1 - (2m.M / r))
          dr = -1 / dt
          dθ = r^2
          dϕ = r^2 * sin(θ)^2
          dtdϕ = zero(r)
          return SVector(dt, dr, dθ, dϕ, dtdϕ)
      end
      ```
    ])
    ],
    [
      #uncover("6-")[
      - No-$bb(Z)_2$
      - Dilaton-Axion
      - Johannsen(-Psaltis)
      - ...
      ]
      #v(1fr)
      #uncover("5-")[
      Sanity check:
      #v(0.5em)
      #text(size: 15pt)[
      ```julia
      using Symbolics, Latexify
      ds = @variables dt, dr, dθ, dϕ, r, θ, M
      comp = metric_components(Schwarzschild(M), (r, θ))
      sum(ds[i]^2 * comp[i] for i in 1:4) |> latexify
      ```
      #v(0.5em)
      $ r^2 dif theta^2 + dif t^2 ( -1 + (2 M)/r ) +  (- dif r^2) / (-1 + (2 M)/r) + sin^2 ( theta ) r^2 dif phi^2 $
      #v(4em)
      ]
    ]
    ]
  )
]

#slide(title: "GPWG")[
  Gradus.jl reproduces tests from *Gold et al., 2020*:
  // figure of the radiative transfer and deflection angle
  #align(center)[
    #image("./figs/gold-2022.png", width: 70%)
  ]
  #grid(
    columns: (65%, 1fr),
    [
      *Mesh geometry* works, but _untested_ with GRMHD (did not have a use case)

      Studying *horizonless compact objects*:
      - Photon rings / lensing bands
      - Pattern speed (Conroy et al., 2023)

      Emergent features:
      - Components added work for *all\* implemented spacetimes*; define new observable, test it on everything
    ],
    [
      #image("./figs/teapot.png", width: 100%)
    ]
  )
]

#slide(title: "Thank you")[
  #v(1em)
  #align(center)[
  #cbox(fill: PRIMARY_COLOR, width: 90%, text(fill: SECONDARY_COLOR)[
    == Summary
    #align(left)[
    *Gradus.jl*: accelerate exploration of new spacetimes
    - *Component modelling* framework
    - *Tested and validated* against a wealth of results in the literature
    - Main applications currently in X-ray astronomy, but easy to adapt
    - Use as a *"try before you buy"* (into developing new software)
    ]
  ])
  ]
  #v(1em)
  #set text(size: 20pt)
  Gradus.jl source code (GPL-3.0):
  - https://github.com/astro-group-bristol/Gradus.jl
  #v(1em)
  Documentation:
  - https://astro-group-bristol.github.io/Gradus.jl/
  #v(1em)
  Source for slides and figures:
  - https://github.com/fjebaker/eht-gpwg-talk-2024
]
