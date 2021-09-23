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
        nw = [(10, 3, 7.66), (20, 4, 9.01), (30, 4,9.63), (40, 5,10.36), (50, 8,12.32), (60, 8,12.63), (70, 8,12.55), (80, 10,14.68), (90, 13,16.76), (100, 13,17.67),
            (110, 15, 18.94), (120, 15, 20.47), (130, 13, 18.23), (140, 17, 21.57), (150, 18, 22.41), (160, 20, 26.23), (170, 19, 24.15), (180, 24, 32.04), (190, 24, 30.27), (200, 25, 32.10)]
        ns = [x[0] for x in nw]
        FS = 8
        with DataPlt(filename="fig1.%s"%tp, figsize=(10,6)) as dp:
            ax = plt.subplot(221)   # treewidth
            cornertex("(a)", ax, offset=(0,0))
            plt.plot(ns, [x[1] for x in nw], label="space complexity (tree width)")
            plt.plot(ns, [x[2] for x in nw], label="time complexity")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("log2 value")

            ax=plt.subplot(222)   # size
            cornertex("(b)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([("counting_sum", r"$\mathbb{R}$"), ("size_max", "T"), ("counting_max", "P1"), ("counting_max2", "P2")]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s (%s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.legend(loc="upper left", fontsize=FS, ncol=2)
            plt.ylim(1e-3,3e3)

            ax = plt.subplot(223)   # IDP
            cornertex("(c)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([("counting_all_(fft)", r"$\mathbb{C}$"), ("counting_all_(finitefield)", r"GF$(p)$"), ("counting_all", "PN")]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s (%s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.legend(loc="upper left", fontsize=FS)
            plt.ylim(1e-3,3e3)

            ax = plt.subplot(224)  # configurations
            cornertex("(d)", ax, offset=(0,0))
            for k,(prefix, l) in enumerate([("config_max", "P1+S1"), ("config_max_(bounded)", "T+bounding"), ("configs_max_(bounded)", "P1+SN+bounding"), ("configs_max", "P1+SN"), ("configs_max2", "P2+SN"), ("configs_all", "PN+SN")]):
                datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s (%s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            plt.yscale("log")
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("time/s")
            plt.legend(loc="upper left", fontsize=FS)
            plt.ylim(1e-3,3e3)
            #plt.plot(twinx(), ns, getindex.(nw, 2), label="tree width",ytickfontsize=12, xticks=:none,yguidefontsize=14,legend=:bottomright, color=:black, ls=:dash, lw=2, ylabel="tree width")

            plt.tight_layout()


    def fig2(self, tp="pdf"):  # independence polynomials
        nw = [(10, 6), (20, 8), (30, 9), (40, 11), (50, 16), (60, 17), (70, 16), (80, 20), (90, 26), (100, 26)]
        ns = [x[0] for x in nw]
        FS = 10
        with DataPlt(filename="fig2.%s"%tp, figsize=(8,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex("(a)", ax, offset=(-0.02,0))
            plt.plot(ns, [x[1] for x in nw])
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("tree width")

            ax=plt.subplot(122)   # size
            cornertex("(b)", ax, offset=(-0.02,0))
            for k,(prefix, l) in enumerate([("counting_sum", r"$\mathbb{R}$"), ("counting_all_(finitefield)", "GF$(p)$"), ("configs_all", "PN+SN")]):
                for device in ["CPU", "GPU"]:
                    datafile = "../../benchmarks/data/maximal-"+prefix+"-r3-"+device+".dat"
                    if os.path.exists(datafile):
                        y = np.loadtxt(datafile) / 1e9
                        plt.plot(ns[:len(y)], y, label="%s (%s)"%(l, device), lw=1.5, color="C%d"%k, ls="--" if device=="GPU" else "-")
            datafile = "../../benchmarks/data/maximal-configs_all-r3bk-CPU.dat"
            y = np.loadtxt(datafile) / 1e9
            plt.plot(ns[:len(y)], y, label="Bron Kerbosch (CPU)", lw=1.5, color="C4", ls="-")
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

fire.Fire(PLT())
