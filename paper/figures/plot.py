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
        nw = [(10, 3), (20, 4), (30, 6), (40, 7), (50, 8), (60, 10), (70, 11), (80, 11), (90, 15), (100, 15),
            (110, 15), (120, 18), (130, 17), (140, 16), (150, 21), (160, 20), (170, 22), (180, 24), (190, 26), (200, 25)]
        ns = [x[0] for x in nw]
        FS = 10
        with DataPlt(filename="fig1.%s"%tp, figsize=(10,6)) as dp:
            ax = plt.subplot(221)   # treewidth
            cornertex("(a)", ax, offset=(0,0))
            plt.plot(ns, [x[1] for x in nw])
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("tree width")

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
            plt.legend(loc="upper left", fontsize=FS)
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
        nw = [(10, 6), (20, 8), (30, 10), (40, 11), (50, 16), (60, 17), (70, 17), (80, 21), (90, 27), (100, 26)]
        ns = [x[0] for x in nw]
        FS = 10
        with DataPlt(filename="fig2.%s"%tp, figsize=(8,3)) as dp:
            ax = plt.subplot(121)   # treewidth
            cornertex("(a)", ax, offset=(0,0))
            plt.plot(ns, [x[1] for x in nw])
            plt.xlabel("number of vertices, $|V|$")
            plt.ylabel("tree width")

            ax=plt.subplot(122)   # size
            cornertex("(b)", ax, offset=(0,0))
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

fire.Fire(PLT())
