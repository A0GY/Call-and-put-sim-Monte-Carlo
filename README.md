# Monte-Carlo European Option Pricer (C++)

> ✏️ Author — **Ali Yassine**   |   BSc Computer Science, University of Reading   |   April 2025

A single-file C++ program that prices at-the-money European **call** and **put** options via
Monte-Carlo simulation, computes a 95 % confidence interval, and benchmarks the result against
the closed-form Black–Scholes price.

---

## 1  Features
| | Detail |
|---|--------|
| **Numerical rigour** | tracks variance → prints *mean ± 1.96·stderr* (95 % CI). |
| **Analytic benchmark** | closed-form Black–Scholes call; reports relative error. |
| **Reproducible** | deterministic seed option (`seed = 42`) or `seed = 0` for true randomness. |
| **Portability** | builds with GCC/Clang/MSVC, C++14+, no external deps. |
| **Ready for parallelism** | drop-in `#pragma omp parallel for` on the loop and compile with `-fopenmp`/`/openmp`. |

---

## 2  Build & run

### GCC / Clang (Linux / WSL / macOS)
```bash
g++ -O3 -std=c++17 option_pricer/main.cpp -o pricer
./pricer            # add -fopenmp for multicore
```

### Visual Studio 2019/2022 (Windows)
1. **Create** a “Console App” project.  
2. **Project  - Properties - C/C++ - Language - C++ Language Standard** → **ISO C++17 (/std:c++17)**.  
3. Add `main.cpp` to the project, **Build** (Ctrl-Shift-B) and **Run** (Ctrl-F5).

---

## 3  Sample output (Release x64, 100 000 paths, seed = 42)
```
The Put fair value is 5.56713  (95% CI ±0.05625)  and expected payoff at maturity is 5.85256
The Call fair value is 10.4741  (95% CI ±0.09624)  and expected payoff at maturity is 11.0111
Black–Scholes analytic call = 10.4506  |  relative error = 0.22484 %
```
*CI width halves when paths quadruple (error ∝ 1/√N).*

---

## 4  Code structure

```
option_pricer/
    main.cpp            # 200 LOC – all logic
README.md
LICENSE                # MIT
```

---

## 5  Numerical background

Under risk-neutral dynamics  

\[
S_T = S_0 \,\exp\!\bigl[\,(r-\tfrac12\sigma^2)T + \sigma\sqrt{T}\,Z\bigr],\quad Z\sim\mathcal N(0,1)
\]

the discounted payoff \( e^{-rT}\max(S_T-K,0)\) is an unbiased estimator of the option’s present value.
