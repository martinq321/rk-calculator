# R_k calculator — large inputs, progress and integer display

This is the exploratory browser calculator, NOT the certified FLINT/Arb hybrid.
The main page has no adjacent-error plot section.

## Update an existing GitHub Pages site

**Replace BOTH `index.html` and `worker.js`** in your existing repository root.
The two `.bin` kernel files have not changed. The updated worker also works with
an older copy of kernels.json from this calculator. To keep the downloadable
Python / Sage / C++ sources and metadata synchronized, upload all this package's
files over the old ones. Leave the Pages settings unchanged: main / root.
After deployment, reload with Ctrl+F5. The worker URL carries a version suffix
to avoid accidentally using the previous 10^12-limited worker.

## New site / local test

Put these files at the repository root (not inside an extra folder):
`index.html`, `worker.js`, `kernels.json`, `kernels_fine.bin`, `kernels_mesh.bin`.
Upload the remaining source and help files as well. Keep `.nojekyll`.
For a local test: `python -m http.server 8000`, then open localhost:8000.
Opening index.html through file:// does not load the worker data correctly.

## Inputs and output

- k = 1,...,6 and 1 < x <= 10^300.
- Enter `10^100`, `10**100`, `1e100`, or a decimal number; commas are accepted.
- The main answer is the numerical approximation rounded to a whole integer.
- Large results also have a compact scientific-notation summary.
- Large results are NOT known to the units digit. Float64 has limited precision;
  formatting an integer, even with BigInt, does not add precision to the calculation.
- The sensitivity diagnostics are NOT error certificates.
- The original plots have not been restored.

## Progress / time bar

The bar advances on fixed kernel bytes downloaded, quadrature nodes processed,
and exact-count sieve stages. Elapsed time is real browser time. No artificial
waiting period is added. A rough remaining download time appears only when enough
transfer measurements exist; no invented countdown is used for instantaneous work.
Fast or cached evaluations can finish before several intermediate states are visible.
Cancel immediately terminates the worker; the next run reloads its data.

## Numerical method

For k>=2, the quadrature approximates
  integral_(1/129)^1 Delta_k(sigma) x^sigma / sigma d sigma.
It is a finite part of the real cut, not a finite Selberg expansion. Tanh-sinh
quadrature, truncated Mobius sums, finite endpoint coverage and float64 errors
remain. The interval (0,1/129) is omitted and is NOT enclosed.

Large-input summation uses exp((sigma-1)*log(x)); the final sum is multiplied by x.
This avoids unnecessarily large intermediate exponentials. The cutoff sensitivity
is now summed from the extra cells directly, so it is not lost by subtracting
nearly identical huge totals. The mesh sensitivity can still miss common errors.

For k=1, the convergent Gram series uses an input-dependent stopping test, not
the old fixed 200-term ceiling. Beyond the saved denominator table, zeta(n+1)-1
is below float64 resolution, so n*zeta(n+1) is replaced by n at working precision.

The new numeric range is an input support limit, NOT a uniform accuracy guarantee.
Full correctly-rounded integers at x=10^100 require a different, multiprecision
and error-certified evaluator; this package does not supply that capability.

## Optional sources

Sage/Python:
  load("rk_calculator.py")
  report = rk_report(k=2, x=10**100, Q=64)
  print(round(report["value"]))

Native C++:
  g++ -O3 -std=c++17 rk_native.cpp -o rk_native
  ./rk_native 2 1e100

The optional SageMathCell page sends code to that external service after the user
chooses to load it. The main browser calculator sends no k or x to a compute server.
The WebAssembly build script remains optional and is not required by index.html.
See TESTING.md for exactly what was executed for this revision.
