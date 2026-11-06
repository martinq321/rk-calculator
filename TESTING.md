# Tests actually performed

- Python numerical core executed (NumPy/SciPy).
- Möbius values checked against SymPy for 1..1000.
- Small exact counts independently checked by factorization.
- Gram-series values compared with 50-digit mpmath Riemann R.
- Quadrature and cutoff comparisons over every plotted point; see validation.json.
- Native C++17 source compiled with g++ and run against the supplied kernels.
- Main calculator UI and Web Worker executed in Chromium using an OFFLINE harness
  serving local data from memory; numerical outputs, input validation, and mobile
  horizontal layout checked. See browser_validation.json. The environment blocked
  live local HTTP navigation, so this is not a deployed HTTP/GitHub test.
- Both vector PDFs rendered and checked.

Not performed: a Sage kernel run; a public SageMathCell evaluation; a live GitHub
Pages deployment; a WebAssembly compilation; certified enclosure of the lower cut.
The claimed testing scope must not be broadened beyond these statements.
