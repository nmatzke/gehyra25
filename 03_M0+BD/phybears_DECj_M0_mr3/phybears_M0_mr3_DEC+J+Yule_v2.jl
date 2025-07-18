#######################################################
# This script shows:
# 1. The calculation of the BioGeoBEARS DEC and DEC+J
#    log-likelihoods, using numeric integration with
#    ClaSSE in Julia. (This requires taking the ClaSSE
#    branch log-likelihood, i.e. lq, and subtracting
#    the birthdeath-process lnL (BD lnL, but without
#    the logfactorial term), as well as the root state
#    frequencies.
#
# 2. Maximum likelihood, in a constrained ClaSSE model,
#    achieving the same ML parameters and lnL as 
#    BioGeoBEARS, for DEC and DEC+J.
#######################################################

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

#######################################################
# BioGeoBEARS results: DEC model
#######################################################
# BioGeoBEARS DEC+J on Gehyra M0_unconstrained ancstates: global optim, 3 areas max. 
# d=3e−04; e=0; j=0.0195; LnL=−261.30

# Basic tree info
numTips = sum(trdf.nodeType .== "tip")
numInternal = sum(trdf.nodeType .== "intern") + sum(trdf.nodeType .== "root")
ttl_tree_length = sum(trdf.brlen[trdf.brlen.>0.0])
birthRate = yuleBirthRate = (numInternal-1) / ttl_tree_length

##############################################
##############################################
# DEC+J ML on Gehyra
##############################################
##############################################

bmo = construct_BioGeoBEARS_model_object()
bmo.init[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "birthRate"] .= birthRate
bmo.est[bmo.rownames .== "deathRate"] .= 0.0
bmo.est[bmo.rownames .== "d"] .= 0.0010
bmo.est[bmo.rownames .== "e"] .= 0.0007
bmo.est[bmo.rownames .== "a"] .= 0.0
bmo.type[bmo.rownames .== "j"] .= "free"
bmo.est[bmo.rownames .== "j"] .= 0.01
numareas = 6
n = numstates_from_numareas(10,6,true)

# CHANGE PARAMETERS BEFORE E INTERPOLATOR
max_range_size = 6 # replaces any background max_range_size=1
inputs = setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
(setup, res, trdf, bmo, files, solver_options, p_Es_v5, Es_tspan) = inputs;
p_Ds_v5 = inputs.p_Ds_v5;





lower = bmo.min[bmo.type .== "free"]
upper = bmo.max[bmo.type .== "free"]
pars = bmo.est[bmo.type .== "free"]
parnames = bmo.rownames[bmo.type .== "free"]
bmo_updater_v1!(inputs.bmo) # works

# Set up DEC ML search
pars = bmo.est[bmo.type .== "free"]
func = x -> func_to_optimize_v7(x, parnames, inputs, p_Ds_v5; returnval="bgb_lnL", printlevel=1)
pars = [0.01, 0.01, 0.01]
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
# BioGeoBEARS DEC+J on Gehyra M0_unconstrained ancstates: global optim, 3 areas max. 
# d=3e−04; e=0; j=0.0195; LnL=−261.30

# NLopt: matches!
# d=0.0003,	e=0.0,	j=0.01972,	Julia_sum_lq=-944.8157, rootstates_lnL=-8.3781,	Julia_total_lnLs1=-953.1938, bgb_lnL=-261.301


# Get the inputs & res:
pars = optx
inputs.bmo.est[inputs.bmo.type .== "free"] .= pars
bmo_updater_v1!(inputs.bmo)
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
resDECj_Yule = deepcopy(inputs.res);
geogfn = lgdata_fn
lnLs_tuple = (total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL, total_loglikelihood=bgb_lnL)
optim_result = build_optim_result(opt, optf, optx, ret)
juliaRes_to_Rdata(inputs.res, trdf, inputs, lnLs_tuple, optim_result, geogfn, trfn; outwd=getwd(), outfns=NaN)
resDECj_Yule_archive = deepcopy((inputs.res, inputs, lnLs_tuple, optim_result, geogfn, trfn));
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
wd = "~/GitHub/gehyra25/03_M0+BD/phybears_DECj_M0_mr3/"  # CHANGE THIS
setwd(wd)
sourceall("/GitHub/PhyBEARS.jl/Rsrc/")
res = PhyBEARS_res_to_BGB_res(outfns=NaN)
resDECj = res  # CHANGE THIS
results_object = res

trfn = res|inputs|trfn				# The pipe characters will be converted to dollar signs
geogfn = res|inputs|geogfn		# The pipe characters will be converted to dollar signs
tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

max_range_size = res|inputs|max_range_size
include_null_range = res|inputs|include_null_range

pdffn = "phyBEARS_Gehyra2_DEC+J+Yule_M0_unconstrained_v1.pdf"  # CHANGE THIS
pdf(pdffn, height=12, width=9)
analysis_titletxt ="PhyBEARS DEC+J+Yule on Gehyra M0_unconstrained"  # CHANGE THIS
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)
"""
rstring = replace(rstring, "|"=>"\$") # replacing the pipes
reval(rstring)

