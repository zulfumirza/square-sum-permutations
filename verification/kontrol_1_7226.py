# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Zülfü Mirza Gür
"""
Gemini'nin "T(k) = O(lambda^k), her lambda > 1,7226" iddiasinin bilgisayar kontrolu.

ONEMLI: Bir program sonsuz bir iddiayi ISPATLAYAMAZ. Bu program sunlari yapar:
  1) Ispatta kullanilan her ara esitsizligi, sonlu bir aralikta tek tek test eder.
  2) Ispattaki "ust toplam" (majorant) gercek T(k)'den buyuk mu, diye bakar.
     Bu, ispatin hangi terimi hangi sinirla degistirdiginin dogru oldugunu gosterir.
  3) Karakteristik denklemlerin koklerini hesaplar (1,8393 / 1,7456 / 1,7226).
Geriye kalan kisim, yani kagit uzerindeki tumevarim, bir cebir adimidir
(dosyanin sonundaki aciklamaya bak).

Calistirma: python3 kontrol_1_7226.py      (K=300 icin yaklasik 30-60 saniye)
"""
import math

K = 300                                    # k ve q icin ust sinir
N = (K * K) // 2 + K + 5

# ------------------------------------------------------------------ diziler
B = [1]                                    # B(0) = 1
for n in range(1, N):
    B.append(sum(B[k * k - n - 1] for k in range(math.isqrt(n) + 1, math.isqrt(2 * n) + 1)))
S = [0]                                    # S(x) = B(0) + ... + B(x-1)
for v in B:
    S.append(S[-1] + v)


def T(q):
    return S[(q * q) // 2] if q >= 0 else 0


def L(q, c):                               # yarim-karenin altindaki c tane B'nin toplami
    return T(q) - S[max(0, (q * q) // 2 - c)]


tum_iyi = True


def kontrol(ad, sonuc):
    global tum_iyi
    tum_iyi &= sonuc
    print(("[TAMAM] " if sonuc else "[HATA ] ") + ad)


# ------------------------------------------------------------------ 1) tam ozdeslik
ok = True
for k in range(3, K):
    e = k % 2
    sag = 1 + sum(T(m) for m in range(2, k) if m % 2 != k % 2) \
        + sum(L(k - 2 * i, i * i - e) for i in range(1, (k - 2) // 2 + 1))
    ok &= (sag == T(k))
kontrol("Tam ozdeslik T(k) = 1 + T(k-1) + T(k-3) + ... + sum L(k-2i, i^2 - e)", ok)

# ------------------------------------------------------------------ 2) son blogu silme (4A)
ok = True
for q in range(4, K):
    x = (q * q) // 2
    for c in range(1, x):
        pmin = math.isqrt(x - c + 1)
        if pmin * pmin < x - c + 1:
            pmin += 1
        if c > 2 * pmin + 1:
            break
        ok &= L(q, c) <= S[((q - 2) ** 2) // 2 + q % 2 + c - 1]
kontrol("4A: L(q,c) <= S(floor((q-2)^2/2) + (q mod 2) + c - 1)", ok)

ok = all(L(k - 2 * i, i * i - k % 2) <= T(k - 2 * i - 1)
         for k in range(6, K) for i in range(1, (k - 2) // 2 + 1)
         if k - 2 * i >= 4 and i * i <= k - 2 * i - 1 + k % 2)
kontrol("4A sonucu: L(k-2i, i^2-e) <= T(k-2i-1)  (i^2 <= k-2i-1+e iken)", ok)

# ------------------------------------------------------------------ 3) 4D: L(q,1) acilimi
ok = True
for q in range(4, K, 2):
    sag = 0
    for r in range(0, q // 2 - 1):
        p0 = q - 2 * r
        n0 = (p0 * p0) // 2 - 1
        for j in range(2, p0):
            if (p0 - j) ** 2 >= n0 + 1:
                sag += B[((p0 - 2 * j) ** 2) // 2 - j * j]
    ok &= (sag == L(q, 1))
kontrol("4D.1: L(q,1) = sum_r sum_{j>=2} B((q-2r-2j)^2/2 - j^2)  (cift q)", ok)

ok1 = all(B[p * p // 2 - j * j] <= T(p) for p in range(4, K, 2) for j in range(1, p) if p * p // 2 - j * j >= 0)
ok2 = all(B[p * p // 2 - j * j] <= T(p - 1) for p in range(4, K, 2) for j in range(2, p)
          if p * p // 2 - j * j >= 0 and j * j <= p - 1)
kontrol("4D.2: B(p^2/2 - j^2) <= T(p)", ok1)
kontrol("4D.2: B(p^2/2 - j^2) <= T(p-1)  (j^2 <= p-1 iken)", ok2)

# ------------------------------------------------------------------ 4) majorant: ispatin ust toplami
# Ispatta T(k)'nin icindeki her terim bir ust sinirla degistiriliyor. Bu ust toplamin
# gercekten T(k)'den buyuk ya da esit oldugunu kontrol ediyoruz (I = sinirsiz).
ok = True
for k in range(6, K):
    e = k % 2
    M = 1 + sum(T(m) for m in range(2, k) if m % 2 != k % 2)       # ana toplam
    for i in range(1, (k - 2) // 2 + 1):
        q = k - 2 * i
        c = i * i - e
        if i == 1 and e == 1:
            continue                                               # tek k'de bu terim 0
        if i == 1 and e == 0:                                      # cift k: L(k-2,1) acilir
            for r in range(0, (k - 2) // 2 - 1):
                p0 = k - 2 - 2 * r
                n0 = (p0 * p0) // 2 - 1
                for j in range(2, p0):
                    if (p0 - j) ** 2 >= n0 + 1:
                        p = p0 - 2 * j
                        M += T(p - 1) if (p >= 4 and j * j <= p - 1) else T(p)
            continue
        if q >= 4 and i * i <= q - 1 + e:
            M += T(q - 1)
        else:
            M += T(q)
    ok &= (T(k) <= M)
kontrol("Ispatin ust toplami (majorant) her k icin T(k)'den buyuk ya da esit", ok)


# ------------------------------------------------------------------ 5) kokler
def kok(f, a, b):
    for _ in range(200):
        m = (a + b) / 2
        if f(a) * f(m) <= 0:
            b = m
        else:
            a = m
    return (a + b) / 2


def Sf(l):
    return (l ** -1 + l ** -5) / (1 - l ** -2)


def G(l):
    return Sf(l) * (Sf(l) + l ** -7 / (1 - l ** -2) ** 2)


tau = kok(lambda l: l ** 3 - l ** 2 - l - 1, 1.5, 2)
lam_c = kok(lambda l: Sf(l) ** 2 - (1 - l ** -4), 1.5, 2)
lam_d = kok(lambda l: G(l) - 1, 1.5, 2)
print(f"\nTribonacci (Teorem 8):        {tau:.10f}")
print(f"Parite baglasimi (4C):        {lam_c:.10f}")
print(f"L(k-2,1) acilimi (4D):        {lam_d:.10f}")
print(f"Altin oran (sani):            {(1 + 5 ** 0.5) / 2:.10f}")
print(f"Gercek oran T({K-1})/T({K-2}): {T(K - 1) / T(K - 2):.10f}")

print("\nSONUC:", "tum kontroller gecti" if tum_iyi else "HATA VAR")

"""
KAGIT UZERINDEKI KISIM (4D tumevarimi), ozet:
Cift indeksler icin T(m) <= C*lambda^m, tek indeksler icin T(m) <= rho*C*lambda^m varsay.
Majoranttaki terimler (yukaridaki 4. kontrol) bu varsayimla toplaninca:
  tek k :  T(k) <= 1 + C*lambda^k * [ S(lambda) + kuyruk ]
  cift k:  T(k) <= 1 + K_I + rho*C*lambda^k * [ S(lambda) + lambda^-7/(1-lambda^-2)^2 + kuyruk ]
burada S(lambda) = (lambda^-1 + lambda^-5)/(1 - lambda^-2). Kuyruklar I buyudukce 0'a gider.
Tumevarimin kapanmasi icin  S/rho < 1  ve  rho*(S + lambda^-7/(1-lambda^-2)^2) < 1  gerekir;
boyle bir rho ancak G(lambda) = S*(S + lambda^-7/(1-lambda^-2)^2) < 1 iken vardir,
yani lambda > 1,7226 iken.
"""
