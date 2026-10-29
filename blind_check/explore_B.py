from series import *
from fractions import Fraction as F
r0=0.3; L=3.6
S=Ser(3,r0,L)
h=F(1,2)
for sgn in ['1','(-1)^m']:
    sg = (lambda m:1) if sgn=='1' else (lambda m:1 if m%2==0 else -1)
    cl=lambda m: S.mono(sg(m),ez=+m,f=[-m,0,0])   # wB^{-m} z^{+m}
    ch=[(1,[0,1,0],lambda m:m,(h,-h)),(-1,[0,0,1],lambda m:-m,(h,-h))]
    I=u1(S,ch,cl,6)
    print('B sign',sgn); show(S,I,['wB','a','b'])
