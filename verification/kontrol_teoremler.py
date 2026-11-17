# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Zülfü Mirza Gür
"""
Kare toplamli blok cozumleri - bilgisayar kontrolleri
Proje: Toplamlari Tam Kare Olan Permutasyonlarin Sayisi Uzerine (OEIS A402123)

Bu program ispatlarin YERINI TUTMAZ. Ispatlarda "sonlu sayida durum
dogrudan kontrol edildi" denen yerleri ve formullerin dogrulugunu kontrol eder.
Calistirma: python3 kontrol_teoremler.py   (birkac saniye surer)
"""
from math import isqrt, exp, sqrt


def blok_sayilari(N, kare_mi):
    """b[n] = 1..n'nin, her blokta en kucuk + en buyuk 'ozel sayi' olacak
    sekilde ardisik bloklara bolunme sayisi. Ozel sayi: kare_mi(s) True."""
    b = [1]  # b[0] = 1: bos dizilis
    for n in range(1, N + 1):
        b.append(sum(b[s - n - 1] for s in range(n + 1, 2 * n + 1) if kare_mi(s)))
    return b


def kare(s):
    return isqrt(s) ** 2 == s


def ucgensel(s):
    m = (isqrt(8 * s + 1) - 1) // 2
    return m * (m + 1) // 2 == s


def tamam(ad, kosul):
    print(("[TAMAM] " if kosul else "[HATA ] ") + ad)
    return kosul


hepsi = True

# ---------------------------------------------------------------- Teorem 3
B = blok_sayilari(25000, kare)
sifirlar = [n for n in range(1, len(B)) if B[n] == 0]
hepsi &= tamam("Teorem 3: B(n)=0 olan n'ler {1,2,4,6,7,9,11,13,17,18,22,26,31,37} (n<=25000)",
               sifirlar == [1, 2, 4, 6, 7, 9, 11, 13, 17, 18, 22, 26, 31, 37])
# taban: 38<=n<=61 icin kalan istisna listesinde olmayan bir kare var
E = set(sifirlar)
taban = all(any(B[s - n - 1] > 0 for s in range(n + 1, 2 * n + 1) if kare(s)) for n in range(38, 62))
hepsi &= tamam("Teorem 3 tabani (38..61)", taban)
# domino: n>=62 icin [n+39, 2n] araliginda kare var
hepsi &= tamam("Teorem 3 domino: n>=62 icin [n+39,2n]'de kare var (n<10^6)",
               all(isqrt(2 * n) ** 2 >= n + 39 for n in range(62, 10 ** 6)))

# OEIS verisi ile karsilastirma
oeis = [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 2, 1, 2, 0, 1, 1, 2, 0, 3, 1, 2, 1, 0, 2, 2,
        3, 2, 4, 0, 3, 1, 2, 1, 5, 2, 6, 3, 3, 3, 4, 1, 2, 7, 4, 7, 3, 9, 3, 6, 2, 4, 3, 9, 2, 11, 6, 13, 6,
        9, 5, 8, 5, 2, 10, 8, 13, 10, 16, 9, 20, 8, 16, 3, 12, 4]
hepsi &= tamam("A402123 DATA satiri programla ayni", B[:len(oeis)] == oeis)

# ---------------------------------------------------------------- Teorem 4
BT = blok_sayilari(5000, ucgensel)
hepsi &= tamam("Teorem 4: ucgensel icin B_T(n)=0 olan n'ler {1,4,10} (n<=5000)",
               [n for n in range(1, len(BT)) if BT[n] == 0] == [1, 4, 10])
hepsi &= tamam("Teorem 4 domino: n>=23 icin n-12 >= 2*sqrt(n)+1 (n<10^6)",
               all(n - 12 >= 2 * sqrt(n) + 1 for n in range(23, 10 ** 6)))

# ---------------------------------------------------------------- Teorem 5
alt = all(B[n] >= exp(0.2 * sqrt(n) - 1.4) - 1e-9 for n in range(38, 89))
hepsi &= tamam("Teorem 5 alt sinir tabani (38..88)", alt)
g89 = 4 * sqrt(2) / (1 + sqrt(1 - 4 * sqrt(2) / sqrt(89)))
hepsi &= tamam("Teorem 5 adim esitsizligi: g(89) < ln2/0.2", g89 < 0.6931471805599453 / 0.2)
hepsi &= tamam("Teorem 5 ust sinir B(n) <= 2^(floor(sqrt(2n))-1) (1<=n<=25000)",
               all(B[n] <= 2 ** (isqrt(2 * n) - 1) for n in range(1, len(B))))

# ---------------------------------------------------------------- Teorem 6
S = [0]
for v in B:
    S.append(S[-1] + v)          # S[x] = B(0)+...+B(x-1)


def T(k):
    return S[(k * k) // 2] if k >= 0 else 0


K = 220  # (K^2)/2 = 24200 <= 25000
ozdeslik = True
for x in range(1, (K * K) // 2):
    top, j = 1, 2
    while (j * j + 1) // 2 <= x - 1:
        top += S[(j * j) // 2] - S[max(0, j * j - x)]
        j += 1
    ozdeslik &= (top == S[x])
hepsi &= tamam("Teorem 6 adim 1: sayma ozdesligi (x < 24200)", ozdeslik)
hepsi &= tamam("Teorem 6 adim 2: indeks esitsizligi (k<220)",
               all((k - i) ** 2 - (k * k) // 2 <= ((k - 2 * i) ** 2) // 2
                   for k in range(3, K) for i in range(1, k - 1) if k - 2 * i >= 0))
F = [0, 1, 1]
while len(F) < K:
    F.append(F[-1] + F[-2])
hepsi &= tamam("Teorem 6: F_k <= T(k) <= 2^(k-2) (2<=k<220)",
               all(F[k] <= T(k) <= 2 ** (k - 2) for k in range(2, K)))

# ---------------------------------------------------------------- Teorem 7 (ucgensel)
BT2 = blok_sayilari(230 * 231 // 4 + 300, ucgensel)
ST = [0]
for v in BT2:
    ST.append(ST[-1] + v)


def tt(k):
    return k * (k + 1) // 2


def TT(k):
    return ST[tt(k) // 2] if k >= 0 else 0


while len(F) < 240:
    F.append(F[-1] + F[-2])
hepsi &= tamam("Teorem 7: F_k <= T_T(k) <= 2^(k-2) (2<=k<230)",
               all(F[k] <= TT(k) <= 2 ** (k - 2) for k in range(2, 230)))
hepsi &= tamam("Teorem 7 i=1 esitligi: t(k-1) - floor(t(k)/2) = floor(t(k-2)/2) (k<5000)",
               all(tt(k - 1) - tt(k) // 2 == tt(k - 2) // 2 for k in range(2, 5000)))

# ---------------------------------------------------------------- Teorem 8
# Tribonacci siniri ve 1,7226 siniri ayri dosyada: kontrol_1_7226.py
hepsi &= tamam("Teorem 8 (ozet): T(k)/T(k-1) <= 1,7226 (k=150..219)",
               all(T(k) / T(k - 1) <= 1.7226 for k in range(150, K)))

print("\nSONUC:", "tum kontroller gecti" if hepsi else "HATA VAR")
