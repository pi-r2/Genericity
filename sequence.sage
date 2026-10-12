#Values
p = 2
m = 4
o = 4
n = 10


Fp = GF(p)


R = PolynomialRing(Fp, ['X'+str(i+1) for i in range(n)], n)
X = R.gens()
print("Parameters:", n,m,o,p)
F = []
for k in range(m):

    if not((n-o+k+1) % 2) :
        fk =    sum([
            X[o+i-k]*X[n-1 - i]   for i in range((n-o+k+1)/2-1)
        ]) + X[0]*X[o-k+(n-o+k+1)/2-1]
    else:
        fk =    sum([
            X[o+i-k]*X[n-1 - i]   for i in range(ceil((n-o+k)/2))
        ])
    F.append(fk)

for f in F :
    print(f)


print(

)
print(jacobian(F,X)[:,::-1])
