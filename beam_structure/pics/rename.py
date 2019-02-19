#!/usr/bin/python

import glob
import os

for t in glob.glob("*.png"):
    on = t.split("_")
    if (len(on) == 5):
        on = on[:-1] + on[4].split('.')
        on[3],on[4] = on[4],on[3]
        nn = on
        nn = nn[:-2] + [nn[4] + '.' + nn[5]]
        nn = '_'.join(nn)
        print t, " -> ", nn
        os.rename(t, nn)
    
