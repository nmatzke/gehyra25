
using Test, PhyBEARS, DataFrames

using Dates									# for e.g. Dates.now(), DateTime
using Distributed						# for e.g. @spawn
using Combinatorics					# for e.g. combinations()
using DataFrames						# for DataFrame()
using DelimitedFiles				# for readdlm()
using NLopt									# seems to be the best gradient-free, box-constrained								
using RCall									# To call R code from Julia

# List each PhyBEARS code file prefix here
using PhyloBits.TrUtils			# for e.g. numstxt_to_df()
using PhyloBits.TreeTable
using PhyBEARS.BGExample
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.Parsers
using PhyBEARS.ModelLikes # e.g. setup_DEC_SSE2
using PhyBEARS.Uppass

"""
# Run with:
cd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/"))
include(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/phybears_M0_mr3_DEC+B+d=e=D_v2.jl"))
"""

setwd(expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/"))

# Input geography
lgdata_fn = expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/geog.data")
geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn)

# Input tree
trfn = expanduser("~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/tree.newick")
tr = readTopology(trfn)
trdf = prt(tr)

#######################################################
# BioGeoBEARS results: DEC model
#######################################################
# BioGeoBEARS DEC on Gehyra M0_unconstrained ancstates: global optim, 3 areas max. 
# d=0.001; e=7e−04; j=0; LnL=−310.93

# Basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
birthRate = yuleBirthRate = (numInternal-1) / ttl_tree_length
d_start = 0.001

bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "e"
bmo.type[bmo.rownames .== "e"] .= "d"
bmo.init[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.init[bmo.rownames .== "deathRate"] .= d_start
bmo.est[bmo.rownames .== "deathRate"] .= d_start
bmo.init[bmo.rownames .== "d"] .= d_start
bmo.est[bmo.rownames .== "d"] .= d_start
bmo.init[bmo.rownames .== "e"] .= d_start
bmo.est[bmo.rownames .== "e"] .= d_start
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
numareas = 6
n = numstates_from_numareas(numareas,numareas,false)
#n = 176            # 10 areas, maxareas 3, 176 states

# Update deathRate based on e (single-area extirpation rate)
# NOTE: This also needs a mu of 0.0 for ranges with more than 1 area
# 
# You can achieve this with an area matrix with areas of size 1.0,
# and a "u" fixed to -100.0. Thus, the e for 1 area is:
#
# 1.0 * 1.0^-100.0 = 1.0
#
# And for 2 areas, it is
# 1.0 * 2.0^-100.0 = 7.888609052210118e-31
bmo.type[bmo.rownames .== "u"] .= "fixed"
bmo.init[bmo.rownames .== "u"] .= -100.0
bmo.min[bmo.rownames .== "u"] .= -100.0
bmo.max[bmo.rownames .== "u"] .= 0.0
bmo.est[bmo.rownames .== "u"] .= -100.0
bmo.type[bmo.rownames .== "u_e"] .= "fixed" # We don't need the extirpation rate to be 0 for large ranges
bmo.init[bmo.rownames .== "u_e"] .= 0.0
bmo.est[bmo.rownames .== "u_e"] .= 0.0
bmo

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
root_age_mult=1.5; max_range_size=6; include_null_range=false; max_range_size=NaN
max_range_size = 6 # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;

numstates = length(inputs.res.likes_at_each_nodeIndex_branchTop[1])
root_age = maximum(trdf[!, :node_age])

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v5, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v5 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, sol_Es_v5=sol_Es_v5);

# Check the interpolator
p_Ds_v5.sol_Es_v5(1.0)
Es_interpolator(1.0)

# Do downpass - slow and fast calculation
(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v5!(res; trdf=trdf, p_Ds_v5=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
# 8.8 seconds, lnLs match BioGeoBEARS
# (8.84, 16, -993.5629081577177, -9.03925286266909, -1002.6021610203868, -310.97299243972896)

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v5, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)
# 0.456 seconds, lnLs match BioGeoBEARS
# (0.456, 16, -993.5629081577177, -9.03925286266909, -1002.6021610203868, -310.97299243972896)




##############################################
##############################################
# DEC+B+d=e=D ML on Gehyra
##############################################
##############################################

bmo = construct_BioGeoBEARS_model_object()
bmo.type[bmo.rownames .== "birthRate"] .= "free"
bmo.type[bmo.rownames .== "deathRate"] .= "birthRate"
bmo.type[bmo.rownames .== "e"] .= "birthRate"
bmo.init[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.init[bmo.rownames .== "deathRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= birthRate
bmo.init[bmo.rownames .== "e"] .= birthRate
bmo.est[bmo.rownames .== "e"] .= birthRate
bmo.est[bmo.rownames .== "d"] .= 0.0010
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.est[bmo.rownames .== "j"] .= 0.0
bmo.type[bmo.rownames .== "j"] .= "fixed"
numareas = 6
n = numstates_from_numareas(numareas,numareas,false)

# Update deathRate based on e (single-area extirpation rate)
# NOTE: This also needs a mu of 0.0 for ranges with more than 1 area
# 
# You can achieve this with an area matrix with areas of size 1.0,
# and a "u" fixed to -100.0. Thus, the e for 1 area is:
#
# 1.0 * 1.0^-100.0 = 1.0
#
# And for 2 areas, it is
# 1.0 * 2.0^-100.0 = 7.888609052210118e-31
bmo.type[bmo.rownames .== "u"] .= "fixed"
bmo.init[bmo.rownames .== "u"] .= -100.0
bmo.min[bmo.rownames .== "u"] .= -100.0
bmo.max[bmo.rownames .== "u"] .= 0.0
bmo.est[bmo.rownames .== "u"] .= -100.0
bmo.type[bmo.rownames .== "u_e"] .= "fixed" # We don't need the extirpation rate to be 0 for large ranges
bmo.init[bmo.rownames .== "u_e"] .= 0.0
bmo.est[bmo.rownames .== "u_e"] .= 0.0
bmo




# CHANGE PARAMETERS BEFORE E INTERPOLATOR
max_range_size = 6 # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=false, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;
p_Ds_v5 = inputs.p_Ds_v5;

lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
#bmo_updater_v1!(inputs.bmo) # works but doesn't update deathRate etc.
bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows) # works

# Set up DEC ML search
pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v5; returnval="lnL", printlevel=1)
pars = [0.01, 0.01]
func(pars)
function func2(pars, dummy_gradient!)
	return func(pars)
end # END function func2(pars, dummy_gradient!)


#######################################################
# Best optimizer found in Julia so far - NLopt
#######################################################
using NLopt
func(pars)
opt = NLopt.Opt(:LN_BOBYQA, length(pars))
ndims(opt)
opt.algorithm
algorithm_name(opt::Opt)
opt.min_objective = func2
opt.lower_bounds = lower::Union{AbstractVector,Real}
opt.upper_bounds = upper::Union{AbstractVector,Real}
opt.lower_bounds
opt.upper_bounds
opt.ftol_abs = 0.001 # tolerance on log-likelihood
(optf,optx,ret) = NLopt.optimize!(opt, pars)
#######################################################
# BioGeoBEARS DEC on Gehyra M0_unconstrained ancstates: global optim, 3 areas max. 
# d=0.001; e=7e−04; j=0; LnL=−310.93

# NLopt: matches!
# d=0.00096,	e=0.00075,	Julia_sum_lq=-993.4815, rootstates_lnL=-9.0747,	Julia_total_lnLs1=-1002.5562, bgb_lnL=-310.9376


# Get the inputs & res:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
#bmo_updater_v1!(inputs.bmo) # works but doesn't update deathRate etc.
bmo_updater_v2(inputs.bmo, inputs.setup.bmo_rows) # works

p_Ds_v5_updater_v1!(p_Es_v5, inputs);
# SEE runtests_ClaSSE_tree_n13_DECj_WORKS.jl
# save_everystep_EQ_false_CAN_MATTER_EVEN_ON_THE_Ds
#inputs.solver_options.save_everystep=false # CAN PRODUCE A -20.9 vs. -20.6 difference!
inputs.solver_options.save_everystep=true	# WORKS!! Can make a difference EVEN ON THE Ds!!

# Solve the Es
print("\nSolving the Es once, for the whole tree timespan...")
prob_Es_v5 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v7_simd_sums, p_Es_v5.uE, Es_tspan, p_Es_v5);
# This solution is an interpolator
sol_Es_v5 = solve(prob_Es_v5, solver_options.solver, save_everystep=true, abstol=solver_options.abstol, reltol=solver_options.reltol);
Es_interpolator = sol_Es_v5;
p_Ds_v7 = (n=p_Es_v5.n, params=p_Es_v5.params, p_indices=p_Es_v5.p_indices, p_TFs=p_Es_v5.p_TFs, uE=p_Es_v5.uE, terms=p_Es_v5.terms, sol_Es_v5=sol_Es_v5);

(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v7!(res; trdf=trdf, p_Ds_v7=p_Ds_v7, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)


#######################################################
# Calculate ancestral range probabilities under the PhyBEARS DEC model,
# ML parameters
#######################################################

# (updates the "res" results object)
uppass_ancstates_v7!(res, trdf, p_Ds_v7, solver_options; use_Cijk_rates_t=false)

# Output the updated res to a series of .Rdata objects that can be read in R
#######################################################
# Ancestral states estimation and plotting
#######################################################
resDEC_BD = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(inputs.res, trdf, inputs, lnLs_tuple, optim_result, geogfn, trfn; outwd=getwd(), outfns=NaN)
resDEC_BD_archive = deepcopy((inputs.res, inputs, lnLs_tuple, optim_result, geogfn, trfn));
"""
# Then run, in R:
"""
# (NOTE: NO COMMENTS are allowed *after* R commands, in the same line)
# $ -- This changes the Julia window to an R window
using RCall

# Load the R string WITHOUT dollar-sign symbols (using pipes instead), then replace
rstring = """
library(ape)
library(cladoRcpp)
library(diversitree)
library(BioGeoBEARS)
wd = "~/GitHub/gehyra25/03_M0+BD/phybears_DEC+B+d=e=D_M0_mr3/"  # CHANGE THIS
setwd(wd)
sourceall("/GitHub/PhyBEARS.jl/Rsrc/")
res = PhyBEARS_res_to_BGB_res(outfns=NaN)
resDEC = res  # CHANGE THIS
results_object = res

trfn = res|inputs|trfn				# The pipe characters will be converted to dollar signs
geogfn = res|inputs|geogfn		# The pipe characters will be converted to dollar signs
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

max_range_size = res|inputs|max_range_size
include_null_range = res|inputs|include_null_range

pdffn = "phyBEARS_Gehyra2_DEC+B+d=e=D_M0_unconstrained_v1.pdf"  # CHANGE THIS
pdf(pdffn, height=12, width=9)
analysis_titletxt ="PhyBEARS DEC+B+d=e=D on Gehyra M0_unconstrained"  # CHANGE THIS
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)
"""
rstring = replace(rstring, "|"=>"\$") # replacing the pipes
reval(rstring)


