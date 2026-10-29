from zfun import *
import time
mp.mp.dps=20
for p in [(0.3,0.25,-2.0,0.8),(0.4,0.6,-1.0,2.0)]:
    t=time.time(); b=ZB_rho(*p); print(p, 'ZB_rho =', mp.nstr(b,15), f'{time.time()-t:.1f}s')
