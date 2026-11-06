# R_k online calculator — exploratory numerical edition

This package evaluates the prescribed real-cut almost-prime functions numerically.
It is NOT the earlier certified FLINT/Arb hybrid and does not provide certificates.
It contains a working static browser calculator, Sage/Python source, a native C++
evaluator using the same prepared quadrature, plots, underlying data and checks.

## Fastest start: GitHub Pages

1. Create a repository, for example `rk-calculator`.
2. Upload the CONTENTS of this directory to its root, keeping the binary kernel files.
3. In Settings > Pages, choose Deploy from a branch, `main`, `/(root)`, then Save.
4. Open the Pages URL GitHub provides. For the account `martinq321` and the repository
   above, the expected URL after deployment is `https://martinq321.github.io/rk-calculator/`.
   This package has not been deployed to that account.

No compiler, Python server or Sage server is needed for the main browser calculator.
It uses JavaScript in a Web Worker and downloads prepared quadrature data. It evaluates
x afresh; it does not interpolate the plotted table. The optional `sage.html` page uses
the external public SageMathCell service, and discloses that service use before loading.

For local testing, run `python -m http.server 8000` in this directory and open
`http://localhost:8000/`. Do not open index.html as file:// because browsers normally
block worker/data loading from local file URLs.

## SageMath / Python

In a local Sage session:

```python
load('rk_calculator.py')
print(rk_report(k=2, x=10**6, Q=64))
export_demo(outdir='rk_demo', xmax=10**6, step=100, Q=128, h=0.05, show=True)
```

For a self-contained public SageMathCell calculation, paste ALL of `rk_sagecell.txt`.
It includes the implementation and sample call; no download/import URL is required.
Change the input call near the bottom. `rk_demo.sage` has the same contents.

In Python, dependencies are numpy, scipy and matplotlib (plots only):

```bash
python -c "from rk_calculator import rk_report; print(rk_report(2,1e6))"
```

The numerical core was executed in Python with NumPy and SciPy. It is compatible with
Sage's Python environment, but was NOT run in a Sage kernel or public SageMathCell here.
No guarantee about the public service's time limits or availability is implied.

## Native C++

From this directory, with a C++17 compiler:

```bash
g++ -O3 -std=c++17 rk_native.cpp -o rk_native
./rk_native 2 1000000
```

On MSYS2/Windows the output can be named `rk_native.exe` and run as `./rk_native.exe`.
The source was compiled and run with g++ in Linux. It reads the supplied little-endian
float64 kernel data; it does not link to FLINT, Arb, GMP or MPFR. This is NOT a port of
the earlier certified executable. It evaluates a prepared quadrature rather than
recomputing the prime-zeta densities in C++.

## Optional C++ -> WebAssembly

Emscripten can compile this small C++ evaluator into WebAssembly:

```bash
bash build_wasm.sh
python -m http.server 8000
```

Then open `wasm_demo.html`. The build emits rk_wasm.js, rk_wasm.wasm and rk_wasm.data.
This optional WebAssembly build was NOT executed here (Emscripten was not installed).
The native C++ build was executed. The main index.html works without this optional build.

GitHub Pages hosts static files; it does not run a native .exe on demand. GitHub Actions
can perform a build, but that is not a persistent per-visitor computation server.
A browser port of the previous FLINT/Arb backend would separately need its dependencies
and parallelism/resource model ported and validated.

## Mathematical target and implementation

Let F_k(s)=[z^k] exp(sum_(m>=1) P(m*s) z^m/m). For real sigma>0, the prescribed
lower-bank value is

P_-(sigma) = sum_(j>=1) mu(j)/j log|zeta(j*sigma)| + i*pi*sum_(j<=1/sigma) mu(j)/j.

The target is R_k(x)=integral_0^1 Delta_k(sigma)*x^sigma/sigma d sigma,
where Delta_k=Im(F_k(sigma-i0))/pi. Prime factors are counted WITH multiplicity.
The recurrence n*f_n=sum_(m=1)^n P_-(m*sigma)*f_(n-m) computes the coefficients.

The default computation for k>=2 uses Q=128 reciprocal cells, so integrates only
1/129 < sigma < 1. Each cell is integrated by tanh-sinh quadrature. The Mobius sum
is truncated after j*m*sigma reaches at least 45; the arithmetic is float64. The
positive-real zeta function is evaluated by scipy.special.zetac, using log1p above
1 to avoid loss of the small tail. The image of a rounded endpoint is never evaluated
at a pole. R_1 is instead computed by its convergent Gram series.

This is NOT a finite Selberg asymptotic expansion, and NOT the incomplete-gamma hybrid.
It directly quadratures a finite part of the real cut. The difference between that
finite part and the full integral is NOT bounded by this implementation.

Changing Q or mesh size gives sensitivity measurements, not certified upper errors.
More precision in the printed output does not supply more rigor. High k and large x
are especially sensitive. The interface's k<=6 and x<=1e12 limits are resource bounds,
not a claim that every value has been validated. Plot checks cover k<=3 and x<=1e6.

## Plots and numerical provenance

`adjacent_data.csv` has the single strictly increasing x-grid used for every curve,
with exact N1,N2,N3 computed by a multiplicity-counting sieve and R1,R2,R3 from the
above numerical methods. No curves are extracted from pictures or generated by AI
image tools; no smoothing or fitted scale is used.

For k=1,2, blue is s_k(x)D_k(x); red is -D_(k+1)(x), D_k=R_k-N_k.
Write t=beta-log(log x), beta=-1.0986500372612+0.24119591642563i and
c2=-0.95484578810794-0.23898284814566i (the manuscript's illustrative constants).
Then s_1=-Re(t) and s_2=-Re((t^2+c2)*conj(t))/(2|t|^2).
These first-zero constants are not independently certified here.

There are 10001 unsmoothed grid points from x=10 to x=10^6, with spacing 100 and the
last endpoint appended. The plot uses Q=128, h=0.05. `validation.json` records cutoff
and mesh refinements, small-count cross-checks and the Gram-series cross-check.
No convergence experiment is represented as a bound for the entire omitted lower cut.


## Primary documentation and references

- GitHub Pages: https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages
- Emscripten: https://emscripten.org/docs/getting_started/Tutorial.html
- SageMathCell embedding: https://sagecell.sagemath.org/static/about.html
- Sage load: https://doc.sagemath.org/html/en/reference/repl/sage/repl/load.html
- SciPy zetac: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.zetac.html
- Construction: https://arxiv.org/abs/2610.11024
- Lucia and later references: https://mathoverflow.net/questions/35927/asymptotic-density-of-k-almost-primes
