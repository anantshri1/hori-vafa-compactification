from series import *
r0=0.3; L=3.6
S=Ser(3,r0,L)
for sgn in ['1','(-1)^m']:
    sg = (lambda m:1) if sgn=='1' else (lambda m:1 if m%2==0 else -1)
    cl=lambda m: S.mono(sg(m),ez=-m,f=[-m,0,0])
    ch=[(1,[0,1,0],lambda m:m,(0,1)),(1,[0,0,1],lambda m:m,(0,1))]
    I=u1(S,ch,cl,6)
    print('A sign',sgn); show(S,I,['w','y1','y2'])
