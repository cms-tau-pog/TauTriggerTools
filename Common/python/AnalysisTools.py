import math
import numpy as np
import scipy
import scipy.optimize
import ROOT

from RootObjects import Histogram, Graph, MultiGraph

def KatzLog(passed, total):
    """Returns 1-sigma confidence interval for a ratio of proportions using Katz-log method."""
    if np.count_nonzero(total) != len(total):
        raise RuntimeError("Total can't be zero")
    if np.count_nonzero(passed < 0) != 0 or np.count_nonzero(total < 0) != 0:
        raise RuntimeError("Yields can't be negative")
    if np.count_nonzero(passed > total) != 0:
        raise RuntimeError("Passed can't be bigger than total")
    if passed[0] == 0 and passed[1] == 0:
        return (0, math.inf)
    if passed[0] == total[0] and passed[1] == total[1]:
        y1 = total[0] - 0.5 if total[0] > 0.5 else total[0] * 0.99
        y2 = total[1] - 0.5 if total[1] > 0.5 else total[1] * 0.99
        # in some sources -1 instead of -0.5 is recommended for y2
    else:
        y1 = passed[0] if passed[0] != 0 else 0.5
        y2 = passed[1] if passed[1] != 0 else 0.5
    n1 = total[0]
    n2 = total[1]
    pi1 = y1 / n1
    pi2 = y2 / n2
    theta = pi1 / pi2
    sigma2 = (1 - pi1) / (pi1 * n1) + (1 - pi2) / (pi2 * n2)
    if sigma2 < 0:
        raise RuntimeError("Invalid inputs: passed={}, total={}".format(passed, total))
    sigma = math.sqrt(sigma2)
    return (theta * math.exp(-sigma), theta * math.exp(sigma))

def weighted_eff_confint_freqMC(n_passed, n_failed, n_passed_err, n_failed_err, alpha=1-0.68, n_gen=100000,
                                max_gen_iters=100, min_stat=80000, seed=42, symmetric=True):
    #print("<weighted_eff_confint_freqMC>:")
    #print(" n_passed = %1.2f +/- %1.2f" % (n_passed, n_passed_err))
    #print(" n_failed = %1.2f +/- %1.2f" % (n_failed, n_failed_err))
    assert n_passed >= 0
    assert n_failed >= 0
    assert n_passed_err >= 0
    assert n_failed_err >= 0
    assert alpha > 0 and alpha < 1
    assert n_gen > 0
    assert max_gen_iters > 0
    assert min_stat > 0
    failed_mc = np.empty(0)
    passed_mc = np.empty(0)
    if seed is not None:
        np.random.seed(seed)
    for gen_iter in range(max_gen_iters):
        new_failed_mc = np.random.normal(n_failed, n_failed_err, n_gen)
        new_passed_mc = np.random.normal(n_passed, n_passed_err, n_gen)
        sel = (new_failed_mc >= 0) & (new_passed_mc >= 0)
        failed_mc = np.append(failed_mc, new_failed_mc[sel])
        passed_mc = np.append(passed_mc, new_passed_mc[sel])
        n_samples = len(passed_mc)
        if n_samples >= min_stat:
            eff = passed_mc / (passed_mc + failed_mc)
            break
    if n_samples < min_stat:
        raise RuntimeError("Unable to estimate confinterval please, increase MC statistics.")
    eff_exp = n_passed / float(n_passed + n_failed)

    if symmetric:
        def coverage(delta_eff):
            x = np.count_nonzero((eff > eff_exp - delta_eff) & (eff < eff_exp + delta_eff))
            return x / float(n_samples)
        opt = scipy.optimize.root_scalar(lambda x: coverage(x) - 1 + alpha, bracket=(0, 1), method='bisect')
        if not opt.converged:
            raise RuntimeError("weighted_eff_confint_freqMC: unable to find a symmetric conf interval.")
        q_down = max(0., eff_exp - opt.root)
        q_up = min(1., eff_exp + opt.root)
        return q_down, q_up
    else:
        eff_up = eff[eff > eff_exp]
        eff_down = eff[eff <= eff_exp]
        frac_up = len(eff_up) / float(n_samples)
        frac_down = len(eff_down) / float(n_samples)
        assert frac_up > 0
        assert frac_down > 0

        def L(alpha_up, return_interal=False):
            alpha_up = min(alpha, max(0, alpha_up))
            alpha_down = (alpha - alpha_up)
            alpha_up_scaled = alpha_up / frac_up
            alpha_down_scaled = min(1., alpha_down / frac_down)

            q_up = np.quantile(eff_up, 1 - alpha_up_scaled) if alpha_up != 0 else 1.
            q_down = np.quantile(eff_down, alpha_down_scaled) if alpha_down != 0 else 0.
            l = q_up - q_down
            if return_interal:
                return l, q_down, q_up
            return l

        opt = scipy.optimize.minimize_scalar(L, bounds=(0, min(alpha, frac_up)), method='Bounded')
        if not opt.success:
            raise RuntimeError("weighted_eff_confint_freqMC: unable to find a conf interval with the minimal size.")

        _, q_down, q_up = L(opt.x, True)
        return q_down, q_up

def ListToStdVector(l, elem_type='string'):
    v = ROOT.std.vector(elem_type)()
    for x in l:
        if elem_type in ['Int_t', 'UInt_t']:
            x = int(x)
        v.push_back(x)
    return v

def RemoveOverflowBins(hist):
    for bin in [ 0, hist.GetNbinsX() + 1 ]:
        hist.SetBinContent(bin, 0)
        hist.SetBinError(bin, 0)

def FixNegativeBins(hist, fix_integral=False, max_rel_shift=0.65):
    has_fixes = False
    integral = hist.Integral()
    if integral <= 0:
        raise RuntimeError("Unable to fix negative bins if integral <= 0.")
    for n in range(hist.GetNbinsX() + 2):
        x = hist.GetBinContent(n)
        if x < 0:
            x_err = hist.GetBinError(n)
            if x + 3.*x_err < 0:
                raise RuntimeError("Yield in bin {} is {} +- {}. Negative bin for which the yield is not statistically"
                                   " compatible with 0 can't be fixed.".format(n, x, x_err))
            hist.SetBinError(n, math.sqrt(x_err ** 2 + x ** 2))
            hist.SetBinContent(n, 0)
            has_fixes = True
    if has_fixes:
        new_integral = hist.Integral()
        total_rel_shift = abs(new_integral - integral) / integral
        if total_rel_shift > max_rel_shift:
            print("total_rel_shift: %1.2f, max_rel_shift %1.2f" % (total_rel_shift, max_rel_shift))
            raise RuntimeError("The overal shift to the integral due to negative bins = {} is above the allowed limit" \
                               " = {}.".format(total_rel_shift, max_rel_shift))
        if fix_integral:
            sf = integral / new_integral
            hist.Scale(sf)

def FixEfficiencyBins(hist_passed, hist_total, remove_overflow=True):
    if remove_overflow:
        RemoveOverflowBins(hist_passed)
        RemoveOverflowBins(hist_total)
    FixNegativeBins(hist_passed)
    FixNegativeBins(hist_total)
    for i in range(hist_total.GetNbinsX() + 2):
        if hist_passed.GetBinLowEdge(i) != hist_total.GetBinLowEdge(i):
            raise ValueError("Histograms passed as function arguments have incompatible binning !!")
        delta = hist_passed.GetBinContent(i) - hist_total.GetBinContent(i)
        if delta > 0:
            if delta > hist_passed.GetBinError(i):
                print("Warning: The number of passed events = {} +/- {} is above the total number events" \
                                   " = {} +/- {} in bin {} [{}, {})." \
                                   .format(hist_passed.GetBinContent(i), hist_passed.GetBinError(i),
                                           hist_total.GetBinContent(i), hist_total.GetBinError(i), i,
                                           hist_total.GetBinLowEdge(i),
                                           hist_total.GetBinLowEdge(i) + hist_total.GetBinWidth(i)))
                print("         Setting bin-content of 'pass' histogram for bin #{} to {}.".format(i, hist_total.GetBinContent(i)))
            hist_passed.SetBinError(i, math.sqrt(hist_passed.GetBinError(i) ** 2 + delta ** 2))
            hist_passed.SetBinContent(i, hist_total.GetBinContent(i))

def dumpHistogram(histName, n_bins, hist_binEdges, hist_binContents, hist_binErrors2):
    if len(hist_binEdges) != (n_bins + 1) or len(hist_binContents) != n_bins or len(hist_binErrors2) != n_bins:
       raise ValueError("Internal error !!")
    print("histogram = %s" % histName)
    print(" bin-contents = ", hist_binContents)
    print(" bin-errors = ", [ math.sqrt(hist_binError2) for hist_binError2 in hist_binErrors2])
    ##print(" bin-error/bin-content = ", [ math.sqrt(hist_binErrors2[i])/hist_binContents[i] if hist_binContents[i] > 0. else 0.5 for i in range(n_bins) ])

def AutoRebinAndEfficiency(hist_passed_a, hist_total_a, hist_passed_b, hist_total_b, max_binError_div_binContent = 0.50):

    ##print("<AutoRebinAndEfficiency>:")

    n_bins = hist_passed_a.GetNbinsX()
    if hist_total_a.GetNbinsX() != n_bins or hist_passed_b.GetNbinsX() != n_bins or hist_total_b.GetNbinsX() != n_bins:
        raise ValueError("Histograms passed as function arguments have incompatible binning !!")

    # CV: convert histograms from ROOT's TH1 to Konstantin's Histogram type
    #    (defined in TauTriggerTools/Common/python/RootObjects.py)
    myhists = [ hist_passed_a, hist_total_a, hist_passed_b, hist_total_b ]   
    for i in range(len(myhists)):
        if type(myhists[i]) != Histogram:
            myhists[i] = Histogram(myhists[i])
    myhist_passed_a = myhists[0]
    myhist_total_a  = myhists[1]
    myhist_passed_b = myhists[2]
    myhist_total_b  = myhists[3]

    ##print("BEFORE rebinning:")
    ##print("#bins = %i" % n_bins)
    ##print(" bin-edges = ", myhist_total_a.edges)
    ##dumpHistogram("data, passed", n_bins, myhist_total_a.edges, myhist_passed_a.values, myhist_passed_a.errors)
    ##dumpHistogram("data, total", n_bins, myhist_total_a.edges, myhist_total_a.values, myhist_total_a.errors)
    ##print("efficiency (data) = ", [ myhist_passed_a.values[i]/myhist_total_a.values[i] for i in range(n_bins) ])
    ##dumpHistogram("mc, passed", n_bins, myhist_total_a.edges, myhist_passed_b.values, myhist_passed_b.errors)
    ##dumpHistogram("mc, total", n_bins, myhist_total_a.edges, myhist_total_b.values, myhist_total_b.errors)
    ##print("efficiency (mc) = ", [ myhist_passed_b.values[i]/myhist_total_b.values[i] for i in range(n_bins) ])

    hists_rebinned_binContents = [ [], [], [], [] ]
    hists_rebinned_binErrors2  = [ [], [], [], [] ]    
    hist_rebinned_binEdges = []
    n_bins_rebinned = 0

    is_unmerged_bin = False

    # merge bins of the original histogram from left to right
    # until sufficient event statistics is accumulated in each bin of each rebinned histogram
    for idx_bin in range(n_bins):

        binEdge = myhist_passed_a.edges[idx_bin]
        if abs(myhist_total_a.edges[idx_bin] - binEdge) > 1.e-1 or abs(myhist_passed_b.edges[idx_bin] - binEdge) > 1.e-1 or abs(myhist_total_b.edges[idx_bin] - binEdge) > 1.e-1:
            raise ValueError("Histograms passed as function arguments have incompatible binning !!")

        if idx_bin == 0:
            hist_rebinned_binEdges.append(myhist_total_a.edges[idx_bin])

        is_sufficient_stats = True
        for idx_hist in range(len(myhists)):
            binContent = myhists[idx_hist].values[idx_bin]
            if len(hists_rebinned_binContents[idx_hist]) < (n_bins_rebinned + 1):
                hists_rebinned_binContents[idx_hist].append(0.)
            hists_rebinned_binContents[idx_hist][n_bins_rebinned] += binContent

            binError2 = myhists[idx_hist].errors[idx_bin] ** 2
            if len(hists_rebinned_binErrors2[idx_hist]) < (n_bins_rebinned + 1):
                hists_rebinned_binErrors2[idx_hist].append(0.)
            hists_rebinned_binErrors2[idx_hist][n_bins_rebinned] += binError2

            # CV: require that all rebinned histograms have non-negative bin-contents
            if not hists_rebinned_binContents[idx_hist][n_bins_rebinned] >= 0.:
                is_sufficient_stats = False
            # CV: require that all rebinned "total" histograms have positive bin-contents
            if not hists_rebinned_binContents[idx_hist][n_bins_rebinned] > 0.:
                is_sufficient_stats = False

            # CV: check sufficient event statistics condition only for "total" histograms
            if (idx_hist == 1 or idx_hist == 3) and (math.sqrt(binError2)/binContent > max_binError_div_binContent):
                is_sufficient_stats = False
        # CV: require that number of events in "passed" histogram is less than or equal to number of events in "total" histogram
        if hists_rebinned_binContents[0][n_bins_rebinned] >= hists_rebinned_binContents[1][n_bins_rebinned] or \
           hists_rebinned_binContents[2][n_bins_rebinned] >= hists_rebinned_binContents[3][n_bins_rebinned]:
            is_sufficient_stats = False
        if is_sufficient_stats:
            hist_rebinned_binEdges.append(myhist_total_a.edges[idx_bin + 1])
            n_bins_rebinned += 1
            is_unmerged_bin = False
        else:
            is_unmerged_bin = True

    # merge events in last two bins in case last bin does not have sufficient event statistics
    if is_unmerged_bin:
        if n_bins_rebinned >= 1:
            # CV: merge bins of the rebinned histogram from right to left 
            #     until sufficient event statistics is accumulated in each bin 
            keep_merging = True
            while keep_merging and n_bins_rebinned >= 1:
                keep_merging = False
                for idx_hist in range(len(myhists)):
                    hists_rebinned_binContents[idx_hist][n_bins_rebinned - 1] += hists_rebinned_binContents[idx_hist][n_bins_rebinned]
                    hists_rebinned_binContents[idx_hist].pop()
                    if hists_rebinned_binContents[idx_hist][n_bins_rebinned - 1] < 0.:
                        keep_merging = True

                    hists_rebinned_binErrors2[idx_hist][n_bins_rebinned - 1] += hists_rebinned_binErrors2[idx_hist][n_bins_rebinned]
                    hists_rebinned_binErrors2[idx_hist].pop()
                if not is_unmerged_bin:
                    print("Warning: Negative number of events encountered in the rightmost bin.")
                    print("         Merging the two rightmost bins...")
                    n_bins_rebinned -= 1
                is_unmerged_bin = False
        else:
            # CV: always create at least one bin, even if the event statistics in that bin is not sufficient
            n_bins_rebinned = 1

        if len(hist_rebinned_binEdges) < (n_bins_rebinned + 1):
            hist_rebinned_binEdges.append(0.)
        hist_rebinned_binEdges[n_bins_rebinned] = myhist_total_a.edges[n_bins]

    if n_bins_rebinned < 2:
        print("Warning: Using max_binError_div_binContent = %1.2f results in a single bin !!" % max_binError_div_binContent)
        print("         Fit of turn-on curve requires at least two bins.")
        print("         Increasing max_binError_div_binContent parameter to %1.2f and trying again..." % (2*max_binError_div_binContent))
        return AutoRebinAndEfficiency(hist_passed_a, hist_total_a, hist_passed_b, hist_total_b, 2*max_binError_div_binContent)

    if len(hist_rebinned_binEdges) != (n_bins_rebinned + 1):
       raise ValueError("Internal error !!")
    ##print("AFTER rebinning:")
    ##print("#bins = %i" % n_bins_rebinned)
    ##print(" bin-edges = ", hist_rebinned_binEdges)
    ##dumpHistogram("data, passed", n_bins_rebinned, hist_rebinned_binEdges, hists_rebinned_binContents[0], hists_rebinned_binErrors2[0])
    ##dumpHistogram("data, total", n_bins_rebinned, hist_rebinned_binEdges, hists_rebinned_binContents[1], hists_rebinned_binErrors2[1])
    ##print("efficiency (data) = ", [ hists_rebinned_binContents[0][i]/hists_rebinned_binContents[1][i] for i in range(n_bins_rebinned) ])
    ##dumpHistogram("mc, passed", n_bins_rebinned, hist_rebinned_binEdges, hists_rebinned_binContents[2], hists_rebinned_binErrors2[2])
    ##dumpHistogram("mc, total", n_bins_rebinned, hist_rebinned_binEdges, hists_rebinned_binContents[3], hists_rebinned_binErrors2[3])
    ##print("efficiency (mc) = ", [ hists_rebinned_binContents[2][i]/hists_rebinned_binContents[3][i] for i in range(n_bins_rebinned) ])

    # compute efficiency and build graph
    graphs_a = MultiGraph(3, n_bins_rebinned)
    graphs_b = MultiGraph(3, n_bins_rebinned)
    for idx_bin_rebinned in range(n_bins_rebinned):
        binEdge_low = hist_rebinned_binEdges[idx_bin_rebinned]
        binEdge_high = hist_rebinned_binEdges[idx_bin_rebinned + 1]
        binCenter = 0.5*(binEdge_low + binEdge_high)

        for hist in [ "a", "b" ]:
            idx_passed = None
            idx_total = None
            graphs = None
            label = None
            if hist == "a":
                idx_passed = 0
                idx_total = 1
                graphs = graphs_a
                label = "data"
            elif hist == "b":
                idx_passed = 2
                idx_total = 3
                graphs = graphs_b
                label = "mc"
            else:
                continue
       
            binContent_passed = hists_rebinned_binContents[idx_passed][idx_bin_rebinned]
            binError_passed = math.sqrt(max(0., hists_rebinned_binErrors2[idx_passed][idx_bin_rebinned]))
            ##print("%s, passed (bin %i): bin-content = %1.2f +/- %1.2f" % (label, idx_bin_rebinned, binContent_passed, binError_passed))
            binContent_total = hists_rebinned_binContents[idx_total][idx_bin_rebinned]
            binError_total = math.sqrt(max(0., hists_rebinned_binErrors2[idx_total][idx_bin_rebinned]))
            ##print("%s, total (bin %i): bin-content = %1.2f +/- %1.2f" % (label, idx_bin_rebinned, binContent_total, binError_total))
            binContent_failed = max(0., hists_rebinned_binContents[idx_total][idx_bin_rebinned] - hists_rebinned_binContents[idx_passed][idx_bin_rebinned])
            binError_failed = math.sqrt(max(0., hists_rebinned_binErrors2[idx_total][idx_bin_rebinned] - hists_rebinned_binErrors2[idx_passed][idx_bin_rebinned]))
            ##print("%s, failed (bin %i): bin-content = %1.2f +/- %1.2f" % (label, idx_bin_rebinned, binContent_failed, binError_failed))
            eff = binContent_passed / binContent_total
            eff_low, eff_high = weighted_eff_confint_freqMC(binContent_passed, binContent_failed, binError_passed, binError_failed)
            graphs.x[idx_bin_rebinned] = binCenter
            graphs.x_error_low[idx_bin_rebinned] = binCenter - binEdge_low
            graphs.x_error_high[idx_bin_rebinned] = binEdge_high - binCenter
            graphs.y[0, idx_bin_rebinned] = binContent_passed
            graphs.y_error_low[0, idx_bin_rebinned] = math.sqrt(hists_rebinned_binErrors2[idx_passed][idx_bin_rebinned])
            graphs.y_error_high[0, idx_bin_rebinned] = graphs_a.y_error_low[0, idx_bin_rebinned]
            graphs.y[1, idx_bin_rebinned] = binContent_total
            graphs.y_error_low[1, idx_bin_rebinned] = math.sqrt(hists_rebinned_binErrors2[idx_total][idx_bin_rebinned])
            graphs.y_error_high[1, idx_bin_rebinned] = graphs_a.y_error_low[1, idx_bin_rebinned]
            graphs.y[2, idx_bin_rebinned] = eff
            graphs.y_error_low[2, idx_bin_rebinned] = eff - eff_low
            graphs.y_error_high[2, idx_bin_rebinned] = eff_high - eff
    ##print("eff_data = ", [ graphs_a.y[2, i] for i in range(n_bins_rebinned) ])
    ##print("eff_mc = ",   [ graphs_b.y[2, i] for i in range(n_bins_rebinned) ])

    return tuple(graphs_a.ToRootGraphs(n_bins_rebinned)) + tuple(graphs_b.ToRootGraphs(n_bins_rebinned))
