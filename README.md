# R_k calculator — k = 1..30, comma-free integers

Exploratory real-cut calculator. NOT the certified FLINT/Arb hybrid.
No adjacent-error plots are shown on the website.

## Updating an existing GitHub Pages site

Replace all FIVE runtime files at the repository root:

    index.html
    worker.js
    kernels.json
    kernels_fine.bin
    kernels_mesh.bin

Both binary tables have changed. Replacing only index.html and worker.js is
NOT sufficient for ranks above 6. Commit all five together. Pages settings do
not change. Refresh the page after deployment (Ctrl+F5 on Windows).
For all downloadable sources to match, upload every file from the full package.

The worker and data URLs carry version/hash queries to reduce stale-cache issues.
The worker rejects metadata that is not from this revision.

## Input and display

- Integer k from 1 through 30.
- 1 < x <= 10^300.
- Input forms include 10^100, 10**100, 1e100 and ordinary decimals.
- Input commas are accepted, but displayed numbers have NO grouping commas.
- The principal output is the numerical approximation rounded to a whole integer.
- Large outputs also include a compact scientific-notation approximation.
- The progress/time bar and Cancel button are retained.
- Whole-integer formatting does not make the trailing digits accurate.

Example: k=2, x=1000000 displays approximately 209965; the exact count is 210035.

## New numerical tables

The tables have been rebuilt through rank 30, not padded with missing values.
They retain Q=128 reciprocal cells and meshes h=.05 and h=.10, but use longer
endpoint tails (tanh-sinh |u|<=5) and preserve distances to both reciprocal
endpoints separately. The regularized zeta pole is evaluated from a local Taylor
series. This avoids losing essential high powers of logarithms when a node
rounds to an endpoint in ordinary floating-point coordinates.

For k>=2, the target remains the finite-cut integral from 1/129 to 1. The lower
cut is OMITTED. Cutoff and mesh changes are sensitivity diagnostics, not upper
bounds. Higher ranks can involve substantial cancellation and may be unreliable
at small x. The numerical function is not the exact almost-prime count.

For k=1 the convergent Gram series is used. Large-input quadrature sums
exp((sigma-1)*log(x)) in float64 and multiplies by x afterwards.

No full-R_k certificate, correctly-rounded full integer, or uniform accuracy
throughout the accepted input range is asserted. BigInt is used for display
only; it does not increase computation precision.

## Local use and optional sources

Serve this directory with `python -m http.server 8000` and open it through HTTP;
file:// does not support the worker/data loading in this page.

Python / Sage:

    load("rk_calculator.py")
    report = rk_report(k=30, x=10**100, Q=64)
    print(round(report["value"]))

The numerical module uses NumPy and SciPy. It was tested in Python, not inside
a Sage kernel. The optional embedded Sage page uses an external service only
when explicitly activated. The main calculator runs entirely in the browser.

C++ (run from the package folder):

    g++ -O3 -std=c++17 rk_native.cpp -o rk_native
    ./rk_native 30 1e100

The optional WebAssembly script has not been executed for this update and is
not required by the website. See TESTING.md for the executed checks.
