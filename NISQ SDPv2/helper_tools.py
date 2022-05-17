# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 01:31:34 2020

@author: Tobias Haug
"""
################################################################################
# Libraries
################################################################################
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams

import operator
from functools import partial, reduce

import numpy as np

import os
################################################################################
flatten = lambda l: [item for sublist in l for item in sublist]


prod = lambda factors: reduce(operator.mul, factors, 1)


def numberToBase(n, b, n_qubits):
    if n == 0:
        return np.zeros(n_qubits,dtype=int)
    digits = np.zeros(n_qubits,dtype=int)
    counter=0
    while n:
        digits[counter]=int(n % b)
        n //= b
        counter+=1
    return digits[::-1]


def get_pauli_mats():
    sx = np.array([
                    [0, 1],
                    [ 1, 0]
                  ], dtype = np.complex128)
    sy = np.array([
                    [0, -1j],
                    [1j, 0]
                  ], dtype=np.complex128)
    sz = np.array([
                    [1, 0],
                    [0, -1]
                  ], dtype=np.complex128)
    Id = np.array([
                    [1, 0],
                    [ 0, 1]
                  ], dtype=np.complex128)
    return [Id, sx, sy, sz]


def plot1D(data, x, saveto="", dataname="", xlabelstring="", ylabelstring="",
           logx=False, logy=False, marker=True, legend=[], ymin=None, ymax=None,
           xmin=None, xmax=None, custommarkersize=7, fillbetween=[],
           customMarkerStyle=[], custom_color_index=[], customplot1DLinestyle=[],
           customlinewidth=[], custom_legend_order=[], custom_error_y=[], 
           xnbins=5, ynbins=4, errorcapsize=3, custom_ticks_x=[], legendcol=1,
           fsizeLegend=19, fig_size=(6,5)):

    rcParams.update({'figure.autolayout': True})
    rcParams['legend.handlelength'] = 1.25
    rcParams['legend.labelspacing'] = 0.2
    rcParams['legend.handletextpad'] = 0.2
    rcParams['legend.borderpad'] = 0.15
    rcParams['legend.borderaxespad'] = 0.25 
    
    #self constructed from color brewer
    colormap=np.array([(56,108,176), (251,128,114),
                       (51,160,44), (253,191,111), (227,26,28),
                       (178,223,138), (166,206,227), (255,127,0),
                       (202,178,214), (106,61,154), (0,0,0)]) / np.array([255., 255., 255.]) 
    
    elements, n_colors = len(data), len(colormap)
    color_index = np.arange(elements) if len(custom_color_index) == 0 else custom_color_index
    
    if len(custom_error_y) == 0:
        custom_error_data_y = [[] for _ in range(elements)]
    else:
        custom_error_data_y = custom_error_y

    fsize = 18
    fsizeLabel = fsize + 12
    # fsizeLegend = fsize + 1
    
    plt.figure(figsize = fig_size)
    if dataname != "":
        pp = PdfPages(saveto + os.path.sep + dataname + '.pdf')

    ax = plt.gca()
    
    lines_list = []

    if len(customplot1DLinestyle) == 0:
        plot1DLinestyle = flatten([["solid", "dashed", "dotted", "dashdot"] \
                                    for i in range(elements // 4 + 1)])
    else:
        plot1DLinestyle = customplot1DLinestyle
    if marker:
        if len(customMarkerStyle) == 0:
            markerStyle = flatten([["o", "X", "v", "^", "s", ">", "<"] \
                                    for i in range(elements // 4 + 1)])
        else:
            markerStyle = customMarkerStyle
    else:
        markerStyle = ["" for i in range(elements)]

    markersize = custommarkersize
    lwglobal = 3
    tickerWidth = 1.2
    minorLength = 4
    majorLength = 8

    if len(customlinewidth) == 0:
        linewidth = [lwglobal for k in range(elements)]
    else:
        linewidth = customlinewidth
        
    dashdef = []
    for i in range(elements):
        if plot1DLinestyle[i] == "solid":
            dashdef.append([1,1])
        elif plot1DLinestyle[i] == "dotted":
            dashdef.append([0.2,1.7])
            linewidth[i] *= 1.7
        elif plot1DLinestyle[i] == "dashed":
            dashdef.append([3,2])
        elif plot1DLinestyle[i] == "dashdot":
            dashdef.append([5,2.5,1.5,2.5])
        else:
            dashdef.append([1,1])
            
    for i in range(elements):
        if len(custom_error_data_y[i]) > 0:
            l = plt.errorbar(x[i], data[i], yerr=custom_error_data_y[i],
                            color=colormap[color_index[i] % n_colors],
                            linewidth=linewidth[i], linestyle=plot1DLinestyle[i],marker=markerStyle[i], ms=markersize,
                            dash_capstyle = "round", capsize=errorcapsize)
        else:
            l, = plt.plot(x[i], data[i], color=colormap[color_index[i] % n_colors],
                            linewidth=linewidth[i], linestyle=plot1DLinestyle[i],marker=markerStyle[i], ms=markersize,
                            dash_capstyle = "round")
        lines_list.append(l)
    
        if(len(custom_error_data_y[i])==0):
            if(plot1DLinestyle[i] in ["dashed", "dotted", "dashdot"]):
                l.set_dashes(dashdef[i])
            
    for i in range(len(fillbetween)):
        ax.fill_between(fillbetween[i][0], fillbetween[i][1], fillbetween[i][2], 
                        color=colormap[fillbetween[i][3] % n_colors], 
                        alpha=fillbetween[i][4])
        
    if xmin != None and xmax != None:        
        ax.set_xlim([xmin, xmax])
        
    if ymin != None and ymax != None:     
         ax.set_ylim([ymin, ymax])

    if logx:
        ax.set_xscale('log')
    if len(custom_ticks_x) > 0:
        ax.set_xticks(custom_ticks_x)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    if logy:
        ax.set_yscale('log')
        
    if xnbins != None and not logx:
        plt.locator_params(axis = 'x', nbins = xnbins)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    if ynbins != None and not logy:
        plt.locator_params(axis = 'y', nbins = ynbins)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        
    plt.tick_params(axis = 'both', which = 'both', width = tickerWidth)
    plt.tick_params(axis = 'both', which = 'minor', length = minorLength)
    plt.tick_params(axis = 'both', which = 'major', length = majorLength)

    ax.get_xaxis().tick_bottom()  
    ax.get_yaxis().tick_left()  
                
    plt.xlabel(xlabelstring)
    plt.ylabel(ylabelstring)
    
    if len(legend) > 0:
        if len(custom_legend_order) > 0:
            lgd = plt.legend([lines_list[i] for i in custom_legend_order], 
                            legend, fontsize=fsizeLegend, ncol=legendcol,
                            columnspacing=0.5)
        else:
            lgd=plt.legend(legend, fontsize=fsizeLegend, ncol=legendcol,
                            columnspacing=0.5)
    
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fsizeLabel)
    for item in ([] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsize)
    
    if dataname != "":
        pp.savefig(plt.gcf())
        pp.close()


def plot2D(data, x, y, dataname="", saveto="", xlabelstring="", ylabelstring="",
            zmin=None, zmax=None, lim_xy=[], Nplot1D=0, plotx1D=[], ploty1D=[],plot1DLinestyle=None, plot1DColor=None, cbarnbins=5, cbar_label=None,
            xnbins=4, ynbins=4):
    if dataname != "":
        pp = PdfPages(saveto + os.path.sep + dataname + ".pdf")
    
    fsize = 18
    fsizeLabel = fsize + 12

    fig = plt.figure(figsize=(6,5))
    ax = plt.gca()
    cmap = "RdYlBu_r"
    shading = "flat"
    deltax = (x[-1] - x[0]) / len(x)
    deltay = (y[-1] - y[0]) / len(y)
    addnum = int(shading == "flat")
    plotx = np.linspace(x[0] - deltax / 2, x[-1] + deltax / 2, 
                        num = len(x) + addnum)
    ploty = np.linspace(y[0] - deltay / 2, y[-1] + deltay / 2,
                        num = len(y) + addnum)
    # plotx = np.append(x, [x[-1] + deltax]) - deltax / 2
    # ploty = np.append(y, [y[-1] + deltay]) - deltay / 2

    vmin = np.amin(data) if zmin == None else zmin
    vmax = np.amax(data) if zmax == None else zmax
        
    plt.pcolormesh(plotx, ploty, data, vmin=vmin, vmax=vmax,
                    cmap=cmap, linewidth=0, rasterized=True, shading=shading,
                    antialiased=False)
    # plt.axis('equal')
    ax.set_aspect(abs((x[-1] - x[0]) / (y[-1] - y[0])))
    
    if Nplot1D > 0 and plotx1D != None:
        #self constructed from color brewer
        colormap = np.array([(56,108,176), (251,128,114),
                           (51,160,44), (253,191,111), (227,26,28),
                           (178,223,138), (166,206,227), (255,127,0),
                           (202,178,214), (106,61,154), (0,0,0)]) / np.array([255., 255., 255.]) 
        n_colors = len(colormap)
        lwglobal = 2
        
        if plot1DLinestyle == None:
            plot1DLinestyle = ["dashed" for i in range(Nplot1D)]
            
        linewidthMod = np.ones(Nplot1D)
        dashdef = []
        for i in range(0,Nplot1D):
            if plot1DLinestyle[i] == "solid":
                dashdef.append([1,1])
                linewidthMod[i] = 1
            elif plot1DLinestyle[i] == "dotted":
                dashdef.append([0.2,1.7])
                linewidthMod[i] = 1.7
            elif plot1DLinestyle[i] == "dashed":
                dashdef.append([3,2])
                linewidthMod[i] = 1
            else:
                dashdef.append([1,1])
                linewidthMod[i] = 1

        linewidthList = [lwglobal * linewidthMod[i] for i in range(Nplot1D)]
        
        if len(ploty1D) == 0:
            ploty1D = [ploty for i in range(Nplot1D)]

        if plot1DColor == None:
            plot1DColor = ["k" for i in range(Nplot1D)]
        else:
            plot1DColor = [colormap[plot1DColor[i] % n_colors] \
                            for i in range(Nplot1D)]
        for i in range(Nplot1D):
            if plot1DLinestyle[i] == "dashed":
                l, = plt.plot(plotx1D[i], ploty1D[i], color=plot1DColor[i], 
                                linewidth=linewidthList[i], 
                                linestyle=plot1DLinestyle[i],
                                dash_capstyle="round")
            else:
                l, = plt.plot(plotx1D[i], ploty1D[i], color=plot1DColor[i],
                                linewidth=linewidthList[i],
                                linestyle=plot1DLinestyle[i],
                                dash_capstyle="round")
            if plot1DLinestyle[i] in ["dashed", "dotted"]:
                l.set_dashes(dashdef[0])
    
    if len(lim_xy) > 0:
        plt.axis([lim_xy[0], lim_xy[1],lim_xy[2], lim_xy[3]])
    else:
        plt.axis([plotx.min(), plotx.max(),ploty.min(), ploty.max()])

    cbar = plt.colorbar(fraction = 0.046, pad = 0.04)
    cbar.locator = MaxNLocator(nbins = cbarnbins)
    cbar.ax.tick_params(labelsize = fsize)
    if cbar_label != None:
        cbar.set_label(cbar_label, rotation=90, size=fsize+3)
                        #, verticalalignment='baseline'
    cbar.update_ticks()

    plt.locator_params(axis = 'x', nbins = xnbins)
    plt.locator_params(axis = 'y', nbins = ynbins)

    plt.xlabel(xlabelstring)
    plt.ylabel(ylabelstring)
    
    for item in ([ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fsizeLabel)
    for item in ([] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fsize)
        
    pp.savefig(plt.gcf(), bbox_inches='tight')
    pp.close()