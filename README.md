# Square-sum permutations

Code accompanying the paper

> Z. M. Gür, *On permutations π of {1,…,n} with i+π(i) a perfect square* (2026).

Let a(n) be the number of permutations π of {1,…,n} with i+π(i) a perfect square for every i
([OEIS A095986](https://oeis.org/A095986)), and B(n) the number of *block solutions*
([OEIS A402123](https://oeis.org/A402123)).

Main results of the paper:

* **Theorem A.** ln S(x) / √x → √2·ln φ, where S(x) = Σ_{m<x} B(m) and φ is the golden ratio.
* **Theorem B.** ln a(n) = ½·n·ln n + O(n).

## Verification scripts

The scripts do **not** replace the proofs. They check the finite cases used in the proofs
and reproduce every numerical claim of the paper.

| Script | Checks | Needs | Time |
|---|---|---|---|
| `verification/kontrol_teoremler.py` | B(n) data vs. OEIS A402123; existence theorem (zeros of B); first bounds; Fibonacci lower bound; triangular analogue; counting identity | Python 3 | seconds |
| `verification/kontrol_1_7226.py` | exact identity (2.1), intermediate inequalities of the earlier λ > 1.7226 bound | Python 3 | ~1 min |
| `verification/kontrol_teorem11.py` | a(1..23) by exact enumeration; Gauss sum bound; scaling function f; W, K, D for n = 500…2500 (λ_min(W), positivity, double stochasticity); explicit threshold n₀ | Python 3 + numpy | ~30 s |

Run e.g. `python3 verification/kontrol_teorem11.py`. Each script prints `[TAMAM]` (OK) or `[HATA]` (error) per check.

## Use of AI tools

This research was carried out with substantial help from AI systems (Claude, Gemini, Grok).
See the section "Use of AI tools" in the paper.

## License

The code in this repository is licensed under the GNU General Public License v3.0 or later (see `LICENSE`).
Copyright (C) 2026 Zülfü Mirza Gür.

The paper itself is not part of this repository; it is distributed separately via arXiv.
