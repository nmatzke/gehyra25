
#######################################################
# Run the PhyBEARS scripts
#######################################################

# These are the ones with include_null_range=false, and the extinction rate for 1-area ranges set by...

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC_M0_mr3/phybears_M0_mr3_DEC+Yule_v2.jl"))

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/phybears_M0_mr3_DEC+J+Yule_v2.jl"))

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D_M0_mr3/phybears_M0_mr3_DEC+B=D_v2.jl"))

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D_M0_mr3/phybears_M0_mr3_DEC+J_B=D_v2.jl"))

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+BD_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+BD_M0_mr3/phybears_M0_mr3_DEC+BD_v2.jl"))

cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+BD_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+BD_M0_mr3/phybears_M0_mr3_DEC+J_BD_v2.jl"))

# birthRate, deathRate=e, no null range
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+e_M0_mr3/phybears_M0_mr3_DEC+B+e_v2.jl"))

# birthRate, deathRate=e, no null range, +J
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B+e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B+e_M0_mr3/phybears_M0_mr3_DEC+J_B+e_v2.jl"))

# birthRate=deathRate=e (no null range)
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D=e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D=e_M0_mr3/phybears_M0_mr3_DEC+B=D=e_v2.jl"))

# birthRate=deathRate=e (no null range), +J
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D=e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D=e_M0_mr3/phybears_M0_mr3_DEC+J_B=D=e_v2.jl"))

# birthRate=e=deathRate (no null range)
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D=e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B=D=e_M0_mr3/phybears_M0_mr3_DEC+B=D=e_v2.jl"))

# birthRate=e=deathRate (no null range), +J
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D=e_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DECj+B=D=e_M0_mr3/phybears_M0_mr3_DEC+J_B=D=e_v2.jl"))


#######################################################
# R code
#######################################################
dirs = c(
"phybears_DEC_M0_mr3",
"phybears_DECj_M0_mr3",
"phybears_DEC+B=D_M0_mr3",
"phybears_DECj+B=D_M0_mr3",
"phybears_DEC+BD_M0_mr3",
"phybears_DECj+BD_M0_mr3",
"phybears_DEC+B+e_M0_mr3",
"phybears_DECj+B+e_M0_mr3",
"phybears_DEC+B=D=e_M0_mr3",
"phybears_DECj+B=D=e_M0_mr3",
"phybears_DEC+B+d=e=D_M0_mr3",
"phybears_DECj+B+d=e=D_M0_mr3")

include_null = c("yes",
"yes",
"yes",
"yes",
"yes",
"yes",
"no",
"no",
"no",
"no",
"no",
"no")


library(ape)
library(cladoRcpp)
library(diversitree)
library(BioGeoBEARS)
sourceall("/GitHub/PhyBEARS.jl/Rsrc/")

i = 1
tmprows = NULL
for (i in 1:length(dirs))
	{
	wd = paste0("~/GitHub/gehyra25/03_M0+BD/", dirs[i], "/")  # CHANGE THIS
	setwd(wd)
	res = PhyBEARS_res_to_BGB_res(outfns=NaN)
	
	lnLs = res$lnLs
	
	pars_to_get = c("d", "e", "j", "birthRate", "deathRate", "u", "u_e", "u_mu")
	param_estimates = res$outputs@params_table[pars_to_get, "est"]
	k = sum(res$outputs@params_table$type == "free")
	tr = ape::read.tree(res$inputs$trfn)
	n = length(tr$tip.label)
	tmprow = c(dirs[i], lnLs, n, k, include_null[i], param_estimates)
	tmprows = rbind(tmprows, tmprow)
	}

tmpdf = as.data.frame(tmprows, stringsAsFactors=FALSE)
names(tmpdf) = c("Model", names(res$lnLs), "n", "k", "null", pars_to_get)
resdf = dfnums_to_numeric(tmpdf)
row.names(resdf) = NULL
resdf

cat(names(resdf), sep='","')
# 
resdf_names = c("Model","total_calctime_in_sec","iteration_number","Julia_sum_lq","rootstates_lnL","Julia_total_lnLs1","bgb_lnL","total_loglikelihood","n","k","null","d","e","j","birthRate","deathRate","u","u_e","u_mu")

subset_names = c("Model", "Julia_total_lnLs1","n","k","null","d","e","j","birthRate","deathRate","u","u_e","u_mu")

subset_names_new = c("Model", "lnL","n","k","null","d","e","j","birthRate","deathRate","u","u_e","u_mu")

resdf2 = resdf[,subset_names]
names(resdf2) = subset_names_new
resdf2

cls.df(resdf2)
resdf2 = unlist_df(resdf2)
cls.df(resdf2)

# Calculated AICc
calc_AICc <- function(lnL, n, k)
	{
	AICc = (-2*lnL) + (2*k*(n/(n-k-1)))
	}

AICc_vals = rep(NA, times=nrow(resdf2))
for (i in 1:nrow(resdf2))
	{
	AICc_vals[i] = calc_AICc(lnL=resdf2$lnL[i], n=resdf2$n[i], k=resdf2$k[i])
	}

AICc_vals = calc_AICc_vals(LnL_vals=resdf2$lnL, nparam_vals=resdf2$k, samplesize=resdf2$n[1])
AICc_vals
AICc = AICc_vals
resdf3 = cbind(resdf2, AICc)

resdf4 = AkaikeWeights_on_summary_table(restable=resdf3, colname_to_use="AICc", add_to_table=TRUE)
resdf4$AICc_wt = resdf4$AICc_wt * 100
cft(resdf4)


# Write out results
wd = "~/GitHub/gehyra25/03_M0+BD/"
setwd(wd)
write.table(x=cft(resdf4), file="PhyBEARS_results_v1.txt", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)







