# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Zülfü Mirza Gür
"""
Teorem 11 (ln a(n) = (1/2) n ln n + O(n)) - bilgisayar kontrolleri
Proje: Toplamlari Tam Kare Olan Permutasyonlarin Sayisi Uzerine

Bu program ispatin YERINI TUTMAZ. Yaptiklari:
  1) a(1..23) degerlerini tam sayimla hesaplar (OEIS A095986 ile karsilastirma).
  2) Lemma 5(1): kuadratik Gauss toplami |G(a,q)|/q <= 1/sqrt(2) (q <= 300), esitlik q = 4'te.
  3) Lemma 1: f(x) * int_0^1 f(y)/(2 sqrt(x+y)) dy = 1 denkleminin cozumu f'yi hesaplar,
     ispattaki [m_f, M_f] sinirlarinin icinde mi diye bakar.
  4) Lemma 2, 3, 6 ve Ana Teorem: belirli n'ler icin W, K, D matrislerini kurar;
     satir toplami sapmasini, ||K||'yi, lambda_min(W)'yi, D'nin pozitif ve cift stokastik
     oldugunu ve max D * sqrt(n) degerini kontrol eder; D'nin verdigi alt siniri yazar.
  5) Ispattaki sabitlerle n_0 esigini hesaplar.
Gereken: numpy.  Calistirma: python3 kontrol_teorem11.py   (yaklasik 1-2 dakika)
"""
import math
from functools import lru_cache
import numpy as np

hepsi = True


def tamam(ad, kosul):
    global hepsi
    hepsi &= bool(kosul)
    print(("[TAMAM] " if kosul else "[HATA ] ") + ad)


kare = lambda v: math.isqrt(v) ** 2 == v

# ------------------------------------------------------------ 1) a(1..23)
def a_tam(n):
    komsu = [[j for j in range(n) if kare(i + j + 2)] for i in range(n)]

    @lru_cache(None)
    def git(i, maske):
        if i == n:
            return 1
        return sum(git(i + 1, maske | 1 << j) for j in komsu[i] if not maske >> j & 1)
    return git(0, 0)


oeis = [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 2, 4, 3, 2, 5, 15, 21, 66, 37, 51]
tamam("a(1..23) tam sayim OEIS A095986 ile ayni", [a_tam(n) for n in range(1, 24)] == oeis)

# ------------------------------------------------------------ 2) Gauss toplamlari
en_buyuk = 0.0
for q in range(2, 301):
    for a in range(1, q):
        if math.gcd(a, q) == 1:
            G = sum(np.exp(2j * np.pi * a * t * t / q) for t in range(q))
            en_buyuk = max(en_buyuk, abs(G) / q)
G4 = abs(sum(np.exp(2j * np.pi * t * t / 4) for t in range(4))) / 4
tamam(f"Lemma 5(1): |G(a,q)|/q <= 1/sqrt2 (q<=300), en buyuk = {en_buyuk:.6f}",
      en_buyuk <= 1 / math.sqrt(2) + 1e-12)
tamam(f"q = 4'te esitlik: |G(1,4)|/4 = {G4:.6f}", abs(G4 - 1 / math.sqrt(2)) < 1e-12)

# ------------------------------------------------------------ 3) f fonksiyonu
Mf, mf = 2 ** (11 / 4), 2 ** (-15 / 4)
izgara = 4000
y = (np.arange(izgara) + 0.5) / izgara
Kc = 1 / (2 * np.sqrt(y[:, None] + y[None, :])) / izgara
u = np.ones(izgara)
for _ in range(3000):
    u = np.sqrt(u / (Kc @ u))          # simetrik olcekleme (sonumlu sabit nokta)
artik = np.max(np.abs(u * (Kc @ u) - 1))
tamam(f"Lemma 1: f bulundu (artik {artik:.1e}), f in [{u.min():.3f}, {u.max():.3f}] "
      f"ispattaki [{mf:.4f}, {Mf:.3f}] icinde", artik < 1e-10 and mf <= u.min() and u.max() <= Mf)
tamam("Lemma 1: f azalmayan", np.all(np.diff(u) >= -1e-12))


def f(x):
    return np.interp(x, y, u)


# ------------------------------------------------------------ 4) W, K, D
CK = Mf * (4 * Mf ** 3 + Mf / 2)
print("\n   n | satir top. sapma | ||K||   | lam_min(W) | max D*sqrt(n) | max|x| | ln a(n) >= ")
for n in [500, 1000, 1500, 2500]:
    i = np.arange(1, n + 1)
    S = i[:, None] + i[None, :]
    A = (np.round(np.sqrt(S)) ** 2 == S).astype(float)
    fv = f(i / n)
    W = fv[:, None] * fv[None, :] * A / math.sqrt(n)
    K = fv[:, None] * fv[None, :] / (2 * math.sqrt(n) * np.sqrt(S))
    s = W.sum(1)
    sapma = np.abs(s - 1).max()
    normK = np.linalg.eigvalsh(K)[-1]
    lam = np.linalg.eigvalsh(W)[0]
    x = np.linalg.solve(np.diag(s) + W, 1 - s)
    D = W * (1 + x[:, None] + x[None, :])
    pozitif = bool(np.all(D[A > 0] > 0))
    stok = np.abs(D.sum(1) - 1).max() < 1e-10 and np.abs(D.sum(0) - 1).max() < 1e-10
    mD = D.max() * math.sqrt(n)
    # per(A) >= per(D) / prod_i max_j D_ij  ve  per(D) >= n!/n^n
    alt = math.lgamma(n + 1) - n * math.log(n) - np.log(D.max(1)).sum()
    print(f"{n:5d} | {sapma:16.4f} | {normK:.4f} | {lam:10.4f} | {mD:13.3f} | {np.abs(x).max():.3f} | {alt:.0f}")
    tamam(f"n={n}: Lemma 2 (sapma <= 2M_f^2/sqrt n), Lemma 3 (||K|| <= 1+C_K/sqrt n), "
          f"lam_min > -0.7472, D pozitif ve cift stokastik",
          sapma <= 2 * Mf ** 2 / math.sqrt(n) and normK <= 1 + CK / math.sqrt(n)
          and lam > -0.7472 and pozitif and stok)
print(f"   (-1/sqrt2 = {-1 / math.sqrt(2):.4f}; ln a(n) alt siniri D'den, kayan nokta hassasiyetinde)")

# ------------------------------------------------------------ 5) n_0
def h(t):
    return np.where(t > 0, np.exp(-1 / np.maximum(t, 1e-300)), 0.0)


xs = np.linspace(-1, 1, 2_000_001)
dx = xs[1] - xs[0]
chi = h(2 * (1 - np.abs(xs))) / (h(2 * (1 - np.abs(xs))) + h(1 - 2 * (1 - np.abs(xs))))
d3 = np.gradient(np.gradient(np.gradient(chi, dx), dx), dx)
Cchi = 4 * (chi.sum() * dx + np.abs(d3).sum() * dx / (2 * math.pi) ** 3)
Cf = 4 * Mf ** 3
Cker, CR = 22 * Mf ** 2, 15
Cloc = 7 * math.sqrt(3) * Cf * Cchi + 2 * math.sqrt(2) * Mf + 9 * Mf * Cchi


def kuv(L, e):                      # n^e, n = 10^L
    return 10 ** (L * e) if L * e > -300 else 0.0


def hata(L):                        # Lemma 6'daki toplam hata
    eker = Cker * 2 ** (1 / 6) * kuv(L, -1 / 120)
    return ((1 / math.sqrt(2)) * (CK * kuv(L, -0.5) + eker) + eker + Mf ** 2 * CR * kuv(L, -1 / 40)
            + 2 * math.sqrt(2) * Mf * Cloc * kuv(L, -3 / 20) + math.sqrt(2) * Cloc ** 2 * kuv(L, -3 / 10))


lo, hi = 18.0, 5000.0
for _ in range(200):
    m = (lo + hi) / 2
    lo, hi = (lo, m) if hata(m) <= 0.04 else (m, hi)
print(f"\nn_0 hesabi (ispattaki sabitlerle, C_chi = {Cchi:.2f}): spektral kosul n >= 10^{hi:.0f} civarinda saglaniyor.")
print("Diger iki kosul (satir toplami, x <= 1/4) n >= 10^18 civarinda saglaniyor; baskin olan spektral kosul.")

print("\nSONUC:", "tum kontroller gecti" if hepsi else "HATA VAR")
