---
title: "Kiddo: High-performance, configurable k-d trees for low-dimensional search in Rust"
tags:
  - Rust
  - nearest-neighbour search
  - spatial indexing
  - k-d trees
  - scientific computing
authors:
  - name: Scott Donnelly
    affiliation: "1"
affiliations:
  - name: Independent Researcher
    index: 1

date: 12 November 2026
bibliography: paper.bib
---

# Summary

Kiddo is an open-source Rust library for fast, flexible nearest-neighbour search
in low-dimensional data. Its modular architecture facilitates configurability for
optimising performance or memory usage for a given use case, or evaluating new
algorithms and layouts. Performance and ergonomics are guiding principles
of the design.

Like other k-d trees, it partitions space to avoid examining distant points
[@bentley1975]. Kiddo combines this pruning with memory layouts designed for
efficient use of processor caches, vectorised distance calculations, prefetching,
and specialised traversal and result collection. These optimisations target
the cost of each search as well as the number of points examined. Author-run
benchmarks report lower query times than other k-d tree implementations in Rust
as well as implementations in other languages [@spatialbench].

The library supports exact and approximate searches, mutable and immutable
datasets, and several coordinate representations. These choices share a common
tree and query interface, so configurations can be compared without rewriting
the surrounding application. Kiddo is distributed through Cargo under the MIT
or Apache-2.0 licence.

# Statement of need

Spatial searches recur throughout scientific computing. Catalogue matching
requires finding nearby sources; point-cloud processing requires local
neighbourhoods; and particle analyses require selecting neighbours within a
cutoff. Repeating a full scan of a large dataset for each query can make these
operations expensive. Using a spatial index library can reduce the work
while leaving the scientific interpretation of the returned neighbours to the
calling program.

Different workloads impose different requirements: an index built once and
queried many times can benefit from an alternative representation from an index
that accepts insertions and removals. Coordinate precision, memory consumption,
query throughput, and startup costs also matter when an index is reused across
analysis runs.

Existing Rust crates and compiled C, C++, and Fortran libraries cover parts of
this space. Kiddo targets the gap with configurable storage and traversal behind
one query interface, so users can choose the trade-off that fits their data.

Kiddo addresses these needs through a common tree and query interface.
Users can select a suitable representation without implementing their own spatial
index. Its scope is
point-based search in low-dimensional coordinate spaces. Suitability and
performance still depend on the data distribution, metric, and workload.
Its intended audience also includes researchers developing spatial data
structures, who need to vary algorithms and memory layouts independently while
retaining a common application interface and evaluation framework.

# State of the field

K-d trees are an established technique for associative search [@bentley1975],
and Friedman et al. established their use for nearest-neighbour search
[@friedman1977]. Broader treatments of multidimensional data structures provide
additional context [@samet2006foundations]. The table below lists notable
packages encountered during Kiddo's development. It is non-exhaustive.

: Related packages for point-based low-dimensional search.

| Package | Language | Queries | Mutability | Design |
|---|---|---|---|---|
| Kiddo | Rust | kNN, radius, approximate nearest-one, best-n, periodic, and batch | mutable and immutable | pluggable strategies, optional SIMD, and serde/rkyv |
| SciPy `cKDTree` [@scipy2020] [@scipy_kdtree] | Python | kNN, radius | immutable | SciPy integration |
| `nanoflann` [@nanoflann2014] | C++ | kNN, radius | mutable and immutable | header-only, dynamic indices |
| `kd-tree` [@terada_kdtree] | Rust | kNN, radius | immutable | custom points, parallel build |
| `kdtree` [@kdtree_crate] | Rust | kNN | mutable | geospatial, custom distance |
| `nabo` [@nabo_rs] [@elseberg2012nabo] | Rust | kNN | immutable | low-dimensional, no_std |
| `neighbourhood` [@neighbourhood_crate] | Rust | radius | immutable | range-only, tunable brute force |
| KDTREE 2 [@kennel2004kdtree2] | Fortran 95/C++ | kNN, radius | immutable | separate Fortran/C++ versions |
| Pkd-tree [@men2025pkdtree] | C++ | kNN, range, and count | mutable | parallel, cache-efficient |

The table separates ecosystem, query scope, mutability, and design. It is
not a performance ranking. Kiddo is the only package here with independent stem
and leaf strategies, so one query interface can be evaluated across storage,
traversal, and result collection.

FLANN provides the algorithmic lineage for nanoflann [@muja2014flann].
Approximate nearest-neighbour search in fixed dimensions has a well-established
theoretical basis [@arya1998approx]. Recent bounds show that exact search can
degrade sharply as dimension grows [@bazzani2026bounds]. Specialised Fortran
implementations continue to appear, including cosmokdtree, an OpenMP-parallel
k-d tree for computational astrophysics [@monllor2026cosmokdtree].

The SPRK-tree study independently evaluated Kiddo 5.0.3 in a comparison of
spherical range-query implementations [@blasius2026sprk]. In its low-dimensional
benchmarks Kiddo was close to the SPRK-tree and placed second on the
lower-dimensional clustering datasets. These results concern the pre-v6
implementation. The Spatial Bench comparison below evaluates the redesigned
Kiddo 6.3.0.

Spatial Bench, also developed by Kiddo author Scott Donnelly, provides
author-run comparisons with Rust, C++, and Python implementations [@spatialbench].
A frozen comparison uses Kiddo 6.3.0, `nanoflann` at revision `7812aa0`, and Rust
`nabo` 0.6.0 [@nabo_rs] for single-threaded exact single-neighbour searches
over three-dimensional uniform `f32` point clouds, using squared-Euclidean
distance on an AMD Ryzen 5 8500GE.
Across nine tree sizes from 65,536 to 16,777,216 points and 1,000 repeated probes,
`nanoflann`'s recorded median query times were approximately 2.2 to 2.6 times
Kiddo's, and `nabo`'s were 2.5 to 3.0 times Kiddo's. The comparison retains fixed
configurations and selects the latest recorded run for each. Timings include
adapter and result handling against prebuilt indexes. Construction is excluded.
The public records and methodology document the configurations and timing
boundaries.

Kiddo's contribution is a Rust spatial-search library whose internal design also
supports data-structure research and performance tuning. The v6 redesign makes
point storage and stem traversal independent extension points behind a shared
tree and query interface. Retrofitting comparable extension points onto a
fixed-layout crate would require breaking changes to its public types. Kiddo's
existing user base made that transition possible without abandoning established
applications. This lets users compare compatible strategies without maintaining separate
search implementations. Comparative performance still
requires matched semantics, data, and build settings.

# Software design

Version 6 was a substantial rewrite developed over approximately a year to
address limitations of earlier versions. Previously, mutable and immutable trees
and floating-point and fixed-point implementations were separated. The
[v6 release notes](https://github.com/sdd-org/kiddo/blob/408e63820bea80811b21ed20c77ff545842d2bd9/CHANGELOG.md#600---2026-07-27)
describe the transition to one generic tree with independent stem and leaf
strategies. The goal was an easier-to-use, more coherent library with a foundation
for correctness, API stability, performance tuning, and experimentation.

## Architecture and extension points

Kiddo exposes storage and traversal as configuration choices behind one query
interface. Users can select the layout that fits their workload without changing
the surrounding application. Shared query semantics, reference tests, and
benchmark infrastructure support comparisons of individual design choices.

## Storage and traversal

Stem strategies include breadth-first Eytzinger layouts [@khuong2017array],
which are usually fastest for larger static arrays, and blocked Donnelly layouts
with configurable block sizes, prefetching, and SIMD traversal variants. Leaf
strategies determine storage layout and mutability, allowing read-optimised and
updatable trees within the same architecture. Their separation exposes
memory-access and construction-versus-query trade-offs for investigation.
Convenience aliases provide defaults for ordinary use. Compile-time strategy
selection avoids runtime trait-object dispatch.

The trait system is itself part of the contribution. Stem and leaf strategies
are separate traits with associated layout and compatibility constraints, and
the tree accepts only compatible combinations through the type system. This
allows block-based traversal, prefetching, and SIMD kernels to share one
implementation without runtime dispatch. The cost is a larger trait surface,
more configuration combinations to validate, and longer compile times.

## Query model

Earlier APIs exposed separate methods for combinations of query options. The
v6 builder composes query family, ordering, inclusive or exclusive boundaries,
periodic images, and result projection, with supported combinations expressed
through type constraints. A new option can be added as another type-state stage
rather than multiplying methods for every combination. Unsupported combinations
fail at compile time, and the same query description can be reused across
strategies or as a target for correctness tests. Results are written into
sorted, threshold, or buffered collectors, and `within_unsorted_visit` avoids
materialising them when the application consumes them directly. Reusing the
same query model also gives correctness checks a common target.

## Metrics and persistence

Generic distance metrics separate stored coordinates from distance output,
permitting wider arithmetic without widening stored points. This makes numerical
trade-offs explicit, particularly for fixed-point data. Exact searches still
follow the chosen arithmetic's rounding and overflow behaviour. Periodic
boundaries support repeating domains, and Rayon supports construction and batch
queries.
Serde and `rkyv` let applications reuse built indexes. Memory-mapped zero-copy
loading removes deserialisation, but a cold page cache still incurs page faults
and storage I/O. Loading is effectively instant only once pages are cached.

## Correctness and performance assurance

Performance is a maintained property: differential and fuzz tests, brute-force
comparisons, regression tests, and scalar paths guard against silent regressions,
supported by continuous benchmark reporting, profiling hooks, and competitor
matrices. SIMD and prefetch kernels are gated, with scalar fallbacks, and unsafe
code is checked across the feature matrix, MSRV, and wasm builds.

# Research impact statement

By 10 November 2026, Kiddo had accumulated over eight million downloads across
all versions on crates.io, exceeding the totals of the other dedicated Rust
k-d tree crates identified in registry searches [@cratesio_downloads2026].

Kiddo has documented uses in astronomy, information-theoretic estimation,
neuron morphology comparison, computational oceanography, and molecular
configurational entropy estimation. Nessie is a
Rust-based friends-of-friends galaxy-group finder with R and Python interfaces.
@lambert2025nessie explicitly identify the Kiddo crate in their method for
estimating survey completeness by counting galaxies around target positions.
This work documents an earlier use of Kiddo rather than adoption of v6.

Carlson Moses Büth's `infomeasure-rs` provides entropy, mutual-information,
transfer-entropy, and related estimators. Its accompanying preprint explicitly
identifies Kiddo as underpinning nearest-neighbour searches
[@buth2026infomeasure]. The inspected implementation uses Kiddo's Euclidean and
Chebyshev queries in its CPU nearest-neighbour estimator paths. Downstream
integration also informs library development: the v6 migration exposed an
empty-input construction bug, reported with a fix and regression tests in
[pull request #388](https://github.com/sdd-org/kiddo/pull/388), and motivated
query-overhead investigations in [issue #394](https://github.com/sdd-org/kiddo/issues/394).

Chris L. Barnes's `nblast-rs` implements neuron morphology comparison and supplies
a Kiddo backend enabled by default [@barnes_nblast]. Its implementation uses
nearest-neighbour queries to estimate local tangent directions and match points
between neurons, with exact and approximate search options. Its manifest targets
the Kiddo 4.2.1-compatible series, demonstrating an established integration with
an earlier API rather than migration to v6.

A public user report additionally describes two-dimensional nearest-neighbour
searches in computational oceanography after upgrading from Kiddo v5 to v6
[@oceanography2026]. This report supports actual use in that field; it does not
identify a particular ocean-circulation model or publication.

Helmut Carter and Tom Kurtzman's `molentropy` estimates molecular configurational
entropy from molecular-dynamics trajectories using internal coordinates and
nearest-neighbour mutual-information expansions [@carter2026molentropy]. Its
implementation calls Kiddo for neighbour searches and handles periodic torsion
coordinates through wrapped query images. The inspected manifest targets the
Kiddo 5.3.1-compatible series. Carter also proposed periodic queries in
[pull request #297](https://github.com/sdd-org/kiddo/pull/297), providing a public
example of research-driven engagement with the library.

Carter confirms that Kiddo supplied nearest-neighbour distances for the entropy
calculations in the Flexi-GIST solvation-thermodynamics work
[@carter2026flexigist], including subsequent revisions (personal communication
supplied by Scott Donnelly, 9 November 2026). This paper-specific connection is
author-confirmed; the inspected preprint does not explicitly name Kiddo.

Kiddo has also served as an experimental baseline in spatial-index research:
@blasius2026sprk cite and evaluate it in their ESA 2026 study of spherical range
queries. This establishes scholarly evaluation alongside the application
integrations described above.

# AI usage disclosure

Generative AI was not used in Kiddo versions before v5. The fundamental v5 and v6
architecture, data-structure designs, and algorithmic ideas are the author's
own. OpenAI Codex assisted with aspects of the implementation alongside
hand-written code, and with analysis of compiler-generated assembly. All changes,
human-written or AI-assisted, go through comprehensive continuous-integration
testing. Test coverage is checked manually during pull-request review and by
automated coverage tools.

OpenAI Codex (GPT-6.1-Sol, high reasoning effort) assisted in drafting this
manuscript, bibliography, and citation metadata, using repository inspection and
primary-source searches. Scott Donnelly repeatedly reviewed, edited, and validated
the manuscript and retains responsibility for core design decisions and technical
claims.

# Acknowledgements

Kiddo was developed independently, with no funding or institutional support,
and outside the auspices of the author's employer.

I thank [Carlson Büth](https://github.com/cbueth),
[Chris Barnes](https://github.com/clbarnes),
[Francisco Couzo](https://github.com/franciscouzo),
[Jules Bertholet](https://github.com/Jules-Bertholet), and
[Perry Kundert](https://github.com/pjkundert) for their contributions to Kiddo,
including work on earlier versions. Special thanks to Carlson Büth for suggesting
a JOSS submission in the first place, his encouragement, and his offer to review
the manuscript. I thank the wider contributor community for development, testing,
issue reports, and discussions.

# References
