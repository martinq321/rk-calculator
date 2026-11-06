# Validation for large-input / progress revision

## Executed

- JavaScript syntax checks for the worker and page script.
- Chromium UI/worker tests using Playwright. The environment blocks top-level URL
  navigation; the page was rendered from its actual HTML and the unchanged worker
  body ran in a Blob with only a relative-fetch origin shim. Its kernel data came
  from a real local HTTP server. This was not a live GitHub Pages deployment test.
- Successful browser calculations for 10^100, 1e100, 10**100, 1,000,000, 10^300,
  a small-x example and a near-one input, covering ranks 1, 2, 3 and 6.
- Seven malformed or out-of-range inputs rejected.
- Progress event percentages were monotone, and completed jobs displayed 100%.
- Slow HTTP data transfer exercised the progress bar, elapsed timer, cancellation
  and restart. The transfer throttle was test-only, not an artificial UI wait.
- Missing-kernel download was reported; retry succeeded after the file returned.
- 390-pixel mobile viewport had no horizontal overflow.
- No JavaScript page errors during the test suite.
- Updated native C++ source compiled with g++ and ran at 1e100 and 1e300.
- Updated Python module executed, including rank 1 and rank 2 at 1e100 and the
  precomputed-kernel vector evaluation through 1e300.

## Numerical comparisons

Browser finite-cut rank 2 at x=10^100:
    2.4821610977539666e98

An independent 40-decimal calculation of the dominant upper real cut gave
    R2_upper / x = 0.0248216109775449363773354331019467796145
It used direct Mobius/zeta boundary values through n=150 and separated the
logarithmic endpoint singularity before quadrature. The browser comparison is
about 2.13e-13 in relative terms. This is a numerical consistency check, not a
rigorous enclosure of omitted intervals, quadrature or all error sources.

Rank-one values were compared with mpmath's Riemann-R evaluator at 45 decimal
working precision. Browser/Python comparisons agreed within 5e-12 relative in
these tests. The figure is a test tolerance, not a certified uniform error bound.

The integer output at x=10^6, k=2 remains 209,965, with exact count 210,035.

Both kernel .bin files are byte-for-byte unchanged from the previous package.
The new algorithm sums scaled exponentials; no new approximation coefficients
are fitted to the example inputs.

## Not claimed

- No complete lower-cut enclosure or correctly rounded full-Rk certificate.
- No promise of uniform accuracy for every supported k and x.
- No meaningful determination of a huge result's units digit.
- No execution in a Sage kernel or the public SageMathCell service in this turn.
- No WebAssembly compilation in this turn.
- No remote GitHub deployment or live-site update performed here.
