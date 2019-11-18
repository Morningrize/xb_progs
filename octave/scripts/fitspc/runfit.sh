#! /bin/bash

echo "./fitspc.m sims/twoplus.kb.dp.xb sims/r$1/events_*MeV_r$1.kb.dp.xb r9xx_cuts_narrow.kb.dp.xb -B r9xx_native_hybrid_bkg.kb.dp.xb -F -o fits/r$1_tp_v01 -x [40:241] -m.lr 1e-3 -m.M 1.5e5" 

./fitspc.m sims/twoplus.kb.dp.xb sims/r$1/events_*MeV_r$1.kb.dp.xb r9xx_cuts_narrow.kb.dp.xb -B r9xx_native_hybrid_bkg.kb.dp.xb -F -o fits/r$1_tp_v01 -x [40:241] -m.lr 1e-3 -m.M 1.5e5
