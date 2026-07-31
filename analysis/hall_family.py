# Independent enumerator of the Prop F2 family (both branches), decades 0..10,
# theta in {0.8, 0.9}. Conventions per the Hall paper: member = distinct x,
# assigned least identity residual; nonsquare x only; exact rational-exponent
# comparison k^q <= x^p; k = 0 members are perfect squares (x^3 = y^2) and
# excluded by the nonsquare rule automatically.
from math import isqrt

XMAX = 10**11

def divisors(n):
    ds, i = [], 1
    while i*i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n//i: ds.append(n//i)
        i += 1
    return ds

def issquare(n):
    if n < 0: return False
    r = isqrt(n); return r*r == n

# w bound: k >= 0.8 w^3 (minus) / 8 w^3 (plus); need k <= x^theta <= XMAX^0.9
# minus: w <= (5/4 * XMAX^0.9)^(1/3); take generous global bound
WBOUND = int((1.25 * XMAX**0.9) ** (1/3)) + 10

fam = {}  # x -> least k
for w in range(1, WBOUND + 1):
    for t in divisors(6*w*w):
        j = 6*w*w // t
        # plus branch
        x = t*t + 4*w
        if x <= XMAX:
            k = 8*w**3 + j*j
            if k < fam.get(x, 1 << 200): fam[x] = k
        # minus branch
        x = t*t - 4*w
        if 2 <= x <= XMAX:
            k = abs(8*w**3 - j*j)
            if 0 < k < fam.get(x, 1 << 200): fam[x] = k

def member(x, k, p, q):
    return (not issquare(x)) and k**q <= x**p

counts = {}
for x, k in fam.items():
    d = len(str(x)) - 1
    for theta, p, q in (("0.8", 4, 5), ("0.9", 9, 10)):
        if member(x, k, p, q):
            counts[(theta, d)] = counts.get((theta, d), 0) + 1

for theta in ("0.8", "0.9"):
    row = [counts.get((theta, d), 0) for d in range(11)]
    print(f"theta={theta}: decades 0..10:", row)
print("published: decade 9 -> 465 / 1568 ; decade 10 -> 832 / 3283")
