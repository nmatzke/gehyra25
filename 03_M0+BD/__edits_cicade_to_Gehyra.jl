phybears_M0_mr3_DEC+B=D_v2.jl

"""
# Run with:
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/phybears_M0_mr3_DEC+J+Yule_v2.jl"))
"""

setwd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/"))

# Input geography
lgdata_fn = expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/geog.data")
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Input tree
trfn = expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/tree.newick")
tr = readTopology(trfn)
trdf = prt(tr)




numareas = 6
n = numstates_from_numareas(10,6,true)

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
max_range_size = 6 # replaces any background max_range_size=1



numareas = 10
n = numstates_from_numareas(10,3,true)

numareas = 6
n = numstates_from_numareas(10,6,true)




wd = "~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/"  # CHANGE THIS



Cicadidae
Gehyra