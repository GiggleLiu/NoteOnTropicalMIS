#!/usr/bin/env python
#The mpl_toolkits.axes_grid module was deprecated in Matplotlib 2.1 and will be removed two minor releases later. Use mpl_toolkits.axes_grid1 and mpl_toolkits.axisartist, which provide the same functionality instead.
import fire
from plotlib import *
import numpy as np
import json
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import os
#from matplotlib.offsetbox import AnchoredText


class PLT(object):
    def fig1(self, tp="pdf"):  # independence polynomials
        ns =[10*i for i in range(1, 26)]
        treewidth = np.loadtxt("../../benchmarks/data/treewidth.dat")
        FS = 8
        with DataPlt(filename="fig1.%s"%tp, figsize=(10,7)) as dp:
            ax = plt.subplot(221)   # treewidth
            cornertex("(a)", ax, offset=(0,0))
            plt.plot(ns, treewidth[:,1], label="space complexity (tree width)", color="C1")
            plt.axhline(y=27, color="C1", ls="--", label="slicing threshold")
            plt.plot(ns, treewidth[:,0], label="time complexity")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("complexity")
            tws = np.arange(0, 45, 5)
            plt.yticks(tws, ["$2^{%d}$"%tw for tw in tws])
            plt.xticks(range(0, 300, 50))
            plt.xlim(0, 260)
            plt.legend(loc="upper left", fontsize=FS, ncol=1)

            ax=plt.subplot(222)   # size
            cornertex("(b)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([("counting_sum", r"counting ISs ($\mathbb{R}$"), ("size_max", "MIS size (T"), ("counting_max", "counting MISs (P1"), ("counting_max2", r"counting MISs and (MIS-1)s (P2")]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s, %s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.legend(loc="upper left", fontsize=FS, ncol=1)
            plt.xticks(range(0, 300, 50))
            plt.xlim(0, 260)
            plt.ylim(1e-3,3e3)

            ax = plt.subplot(223)   # IDP
            cornertex("(c)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([
                ("counting_all", "direct (PN"),
                ("counting_all_(finitefield)", r"finite field fitting (GF$(p)$"),
                ("counting_all_(fft)", r"fourier transform ($\mathbb{C}$"),
                ]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s, %s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.xticks(range(0, 300, 50))
            plt.xlim(0, 260)
            plt.legend(loc="upper left", fontsize=FS)
            plt.ylim(1e-3,3e3)

            ax = plt.subplot(224)  # configurations
            cornertex("(d)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([("config_max", "one MIS (P1+S1"), ("config_max_(bounded)", "one MIS (T+bounding"), ("configs_max_(bounded)", "MISs (P1+SN+bounding"), ("configs_max", "MISs (P1+SN"), ("configs_max2", "MISs and (MIS-1)s (P2+SN"), ("configs_all", "ISs (PN+SN")]):
                datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s, %s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.xticks(range(0, 300, 50))
            plt.xlim(0, 260)
            plt.legend(loc="upper left", fontsize=FS)
            plt.ylim(1e-3,3e3)
            #plt.plot(twinx(), ns, getindex.(nw, 2), label="tree width",ytickfontsize=12, xticks=:none,yguidefontsize=14,legend=:bottomright, color=:black, ls=:dash, lw=2, ylabel="tree width")

            plt.tight_layout()


    def fig2(self, tp="pdf"):  # independence polynomials
        nw = [(10, 6, 9.49), (20, 8, 12.40), (30, 9, 13.35), (40, 11, 15.14), (50, 16, 21.26), (60, 17, 22.58), (70, 16, 21.50), (80, 20, 27.48), (90, 26, 35.04), (100, 26, 35.06)]
        ns = [x[0] for x in nw]
        FS = 8
        with DataPlt(filename="fig2.%s"%tp, figsize=(8,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex("(a)", ax, offset=(-0.02,0))
            plt.plot(ns, [x[1] for x in nw], label="space complexity (tree width)")
            plt.plot(ns, [x[2] for x in nw], label="time complexity")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("complexity")
            tws = np.arange(0, 36, 5)
            plt.yticks(tws, ["$2^{%d}$"%tw for tw in tws])
            plt.legend(loc="upper left", fontsize=FS, ncol=1)

            ax=plt.subplot(122)   # size
            cornertex("(b)", ax, offset=(-0.02,0))
            for k,(prefix, l) in enumerate([("counting_sum", r"counting of maximal ISs ($\mathbb{R}$"), ("counting_all_(finitefield)", "maximal independence polynomial (GF$(p)$"), ("configs_all", "maximal ISs (PN+SN")]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/maximal-"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s, %s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            datafile = "../../benchmarks/data/maximal-configs_all-r3bk-CPU.dat"
            y = np.loadtxt(datafile) / 1e9
            plt.plot(ns[:len(y)], y, label="maximal ISs (Bron Kerbosch, CPU)", lw=1.5, color="C4", ls="-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.legend(loc="upper left", fontsize=FS)
            plt.ylim(1e-3,3e3)

            plt.tight_layout()

    def fig3(self, tp="pdf"):  # independence polynomials
        with DataPlt(filename="fig3.%s"%tp, figsize=(7,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex2("(a)", ax, offset=(0,0))
            seed = 3
            datafile = "../../democode/data/r3_n200_hamming_%d.txt"%seed
            y = np.loadtxt(datafile)
            n = np.max(np.where(y > 0)[0])
            plt.bar(range(n), y[:n], lw=1.5, color="C4", ls="-")
            plt.xlabel("Hamming distance")
            plt.ylabel("Count")
            plt.xlim(0, 200)
            ax.ticklabel_format(style="scientific", scilimits=(-3, 3))

            ax = plt.subplot(122)   # treewidth
            cornertex2("(b)", ax, offset=(0,0))
            seed = 8
            datafile = "../../democode/data/r3_n200_hamming_%d.txt"%seed
            y = np.loadtxt(datafile)
            n = np.max(np.where(y > 0)[0])
            plt.bar(range(n), y[:n], lw=1.5, color="C4", ls="-")
            plt.xlabel("Hamming distance")
            plt.xlim(0, 200)

            plt.tight_layout(h_pad=30.0)

    def fig4(self, tp="pdf"):  # independence polynomials
        with DataPlt(filename="fig3.%s"%tp, figsize=(7,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex2("(a)", ax, offset=(0,0))
            nv = 100
            seed = 3
            datafile = "../../democode/data/r3_n%d_hamming_%d.txt"%(nv,seed)
            y = np.loadtxt(datafile)
            n = np.max(np.where(y > 0)[0])
            plt.bar(range(n), y[:n], lw=1.5, color="C4", ls="-")
            plt.xlabel("Hamming distance")
            plt.ylabel("Count")
            plt.xlim(0, nv)
            ax.ticklabel_format(style="scientific", scilimits=(-3, 3))

            ax = plt.subplot(122)   # treewidth
            cornertex2("(b)", ax, offset=(0,0))
            seed = 5
            datafile = "../../democode/data/r3_n%d_hamming_%d.txt"%(nv,seed)
            y = np.loadtxt(datafile)
            n = np.max(np.where(y > 0)[0])
            plt.bar(range(n), y[:n], lw=1.5, color="C4", ls="-")
            plt.xlabel("Hamming distance")
            plt.xlim(0, nv)

            plt.tight_layout(h_pad=30.0)

    def fig5(self, tp="pdf"):  # 
        def plot_and_plot(which, label, errorbase):
            filename = "../../project/data/%s/entropyconstant_summary.dat"%which
            mean = np.loadtxt(filename)[:,0]
            error = np.sqrt(np.loadtxt(filename)[:,1]/1000*1000)
            ns = np.arange(2,len(mean)+2)
            if errorbase:
                plt.errorbar(ns, mean, label=label, yerr=error)
            else:
                plt.plot(ns, mean, label=label)
        with DataPlt(filename="fig5.%s"%tp, figsize=(8,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex("(a)", ax, offset=(-0.02,0))
            plot_and_plot("square", "Square lattice", False)
            #plot_and_plot("square-0.8", "Square lattice (0.8 filling)", True)
            plt.xlim(0, 42)
            plt.ylim(1.3, 1.7)
            plt.legend()
            plt.xlabel("lattice size L")
            plt.ylabel(r"$F(L,L)^{1/{L^2}}$")
            ax=plt.subplot(122)   # size
            cornertex("(b)", ax, offset=(-0.02,0))
            plot_and_plot("diag", "King's graph", False)
            plot_and_plot("diag-0.8", "King's graph (0.8 filling)", True)
            plt.ylim(1.3, 1.7)
            plt.xlim(0, 34)
            plt.xlabel("lattice size L")
            plt.ylabel(r"$F(L,L)^{1/{L^2}}$")
            plt.legend()
            plt.tight_layout()

    def preprocess_entropy(self, which, nmax):
        meanvar = np.zeros((nmax-1, 2))
        for n in range(2, nmax+1):
            filename = "../../project/data/%s/entropyconstant-%d.dat"%(which, n)
            data = np.loadtxt(filename) ** (1/n**2)
            meanvar[n-2,0] = np.mean(data)
            meanvar[n-2,1] = np.var(data)
        fout = "../../project/data/%s/entropyconstant_summary.dat"%(which,)
        np.savetxt(fout, meanvar)

fire.Fire(PLT())
