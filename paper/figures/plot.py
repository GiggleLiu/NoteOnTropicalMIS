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
            for k,(prefix, l) in enumerate([("counting_sum", "total size"), ("size_max", "max size"), ("counting_max", "max counting"), ("counting_max2", "max2 counting")]):
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
            for k,(prefix, l) in enumerate([("counting_all_(fft)", "IDP (FFT)"), ("counting_all_(finitefield)", "IDP (finite field)"), ("counting_all", "IDP (polynomial)")]):
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
            for k,(prefix, l) in enumerate([("config_max", "single MIS"), ("config_max_(bounded)", "single MIS (bounding)"), ("configs_max_(bounded)", "all MIS (bounding)"), ("configs_max", "all MIS"), ("configs_max2", "all MIS and MIS-1"), ("configs_all", "all IS")]):
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


fire.Fire(PLT())
