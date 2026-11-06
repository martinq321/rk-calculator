# Checks for ranks30-ungrouped-20261106-v1

## Executed functional checks

- Rebuilt BOTH binary quadrature tables for ranks 1..30.
- JavaScript syntax checks passed for the page and worker.
- Chromium tests: 120 actual-worker calculations, every rank 1..30 at
  x=2, 10^6, 10^100 and 10^300, all produced finite numerical results.
- UI checks covered ranks 2, 7, 15 and 30 and 10^100 / 10^300 inputs.
- Output and exact-count comparison have no comma separators.
- Rank 31 is rejected by form validation.
- Progress events were monotone; successful jobs displayed 100%.
- Cancellation terminated the running worker; a subsequent calculation succeeded.
- A 390-pixel mobile layout showed no horizontal overflow.
- No JavaScript page errors were recorded.
- The C++ source compiled and ran with k=30, x=1e100, agreeing with the browser.
- The Python module and a further h=.025 quadrature were executed.

Browser environment limitation: top-level HTTP navigation is blocked by the
runtime. The tests rendered the unchanged page HTML with set_content and used
the unchanged worker body in a Blob. A test-only fetch wrapper resolved relative
paths to a local HTTP server. This was not a live GitHub Pages deployment test.
A delay in the test-only fetch wrapper allowed cancellation to be exercised;
no artificial delay is included in the shipped calculator.

## Numerical consistency checks

The old |u|<=3.15 endpoint coverage was inadequate for high-rank log powers.
The rebuilt tables retain pole distances even when sigma rounds to an endpoint
and extend |u| to 5. Both meshes share this handling.

Additional h=.05 versus h=.025 checks covered k=2,6,10,15,20,25,30 at
x=2, 10^4, 10^6, 10^10, 10^30, 10^100, 10^300. At k=30 and x=10^100,
the relative change was about 2.34e-9. At k=30 and x=10^6 the absolute
change was about 1.12e-11. These are observed changes, NOT error bounds.
A change smaller than 0.5 does NOT certify integer rounding.

The browser value for k=30, x=10^100 was about 7.676912979e92.
For k=2, x=10^6 the integer display is 209965, and N_2(10^6)=210035.

The regularized log|zeta(1+epsilon)| was compared with 160-decimal mpmath
calculations for epsilon of both signs and magnitudes 0.049, 1e-5, 1e-20 and
1e-80; the converted float64 results agreed in these tests.

Full tables of results: browser_validation.json and numerical_validation.json.

## Not established

- No complete lower-cut enclosure or correctly rounded full-R_k certificate.
- No uniform high-rank accuracy guarantee; cancellation can be important.
- No validated units digit for large displayed integers.
- No independent proof of correctness from numerical refinement alone.
- No execution in Sage or the public SageMathCell service for this revision.
- No WebAssembly build or live GitHub deployment performed in this revision.
