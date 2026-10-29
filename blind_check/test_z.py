from zfun import *
import time
for p in [(0.3,0.25,-2.0,0.8),(0.4,0.6,-1.0,2.0)]:
    t=time.time(); a=Z_CP1(*p); t1=time.time()
    print(p,'Z_CP1 =',mp.nstr(a,15), f'({t1-t:.1f}s)')
