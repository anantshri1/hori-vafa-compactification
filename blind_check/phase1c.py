"""Phase 1(c): exact x-expansions of I_A and I_B (mapped), r kept symbolic in exponents."""
from series import *
from fractions import Fraction as F
r0=0.3; L=4.4
S=Ser(3,r0,L)   # fugacity vars: (w=omega, y1, y2)
clA=lambda m: S.mono(1,ez=-m,f=[-m,0,0])
chA=[(1,[0,1,0],lambda m:m,(0,1)),(1,[0,0,1],lambda m:m,(0,1))]
IA=u1(S,chA,clA,8)
q=F(1,4); h=F(1,2)
clB=lambda m: S.mono(1,ez=+m,f=[0,-m,m])          # (y1/y2)^{-m} z^{+m}
chB=[(1,[h,-q,-q],lambda m:m,(h,-h)),(-1,[h,-q,-q],lambda m:-m,(h,-h))]
IB=u1(S,chB,clB,8)
print("I_A expansion:"); show(S,IA,['w','y1','y2'])
print("exact equality of truncated series:", IA==IB, " #terms",len(IA))
# flux-sector bookkeeping for B: which B gauge flux produces each term
for target in [ (F(1),F(-1),F(0),F(1),F(-1,2),F(-1,2)) ]:
    pass
print("\nper-flux origin (terms up to w=3.5):")
def sector(chirals,cl,m):
    t=cl(m)
    for (zc,f,fl,R) in chirals: t=S.mul(t,S.chir(zc,f,fl(m),R))
    return {k:v for k,v in t.items() if k[2]==0}
for m in range(-3,4):
    a=sector(chA,clA,m); b=sector(chB,clB,m)
    print(f"--- gauge flux m={m}:  A-terms:"); show(S,a,['w','y1','y2'],3.5)
    print(f"    B-terms:"); show(S,b,['w','y1','y2'],3.5)
