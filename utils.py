import os
import sys 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm 
from mpl_toolkits.axes_grid1 import make_axes_locatable


class graphTool():
    def __init__(self, figshape, figsize):
        plt.style.use('ggplot')
        mpl.rcParams['axes.unicode_minus'] = False
        self.fig, self.axes = plt.subplots(*figshape, figsize=figsize, squeeze = False)
    
    def params(func):
        def wrapper(*args, **kwargs):
            label = ''
            c = 'r'
            s = None
            title = ''
            for key, val in kwargs.items():
                if key == 'label': label = val
                if key == 'c': c = val
                if key == 's': s = val
                if key == 'title': title = val
            func(*args, **kwargs)
        return wrapper 
        
    @params
    def plot(self, ax_num, x, y = None, label='', c='r', title= ''):
        self.axes[ax_num].plot(x, y, label = label, c = c)
        self.axes[ax_num].legend()
        if title != '':
            self.axes[ax_num].set_title(title)
    @params
    def scatter(self, ax_num, x, y, label='', c='r', s = None, title=''):
        self.axes[ax_num].scatter(x, y, label = label, c = c, s = s)
        self.axes[ax_num].legend()
        if title != None:
            self.axes[ax_num].set_title(title)
    
    def imshow(self, ax_num, image, grid = False, **kwargs):
        cmap = cm.bwr
        no_cbar_ticks = False
        vmax = np.max(image)
        vmin = np.min(image)
        title = None
        for key, val in kwargs.items():
            if key == 'cmap': cmap = val
            if key == 'no_cbar_ticks': no_cbar_ticks = val
            if key == 'vmax': vmax = val
            if key == 'vmin': vmin = val
            if key == 'title': title = val

        im = self.axes[ax_num].imshow(image, vmax=vmax, vmin=vmin, cmap=cmap)
        self.axes[ax_num].set_title(title)
        self.axes[ax_num].grid(False)

        divider = make_axes_locatable(self.axes[ax_num])
        cax = divider.append_axes('left', size='5%', pad=0.1)
        self.axes[ax_num].yaxis.set_ticks_position('right')
        cbar = self.fig.colorbar(im, cax=cax, orientation='vertical')
        cax.yaxis.set_ticks_position('left')
        if no_cbar_ticks == True: cbar.set_ticks([])

    def save(self, path):
        plt.savefig(path)
        
    def show(self):
        plt.show()
        plt.close('all')


def read_excel(path, encoding = 'utf-8', header=[1,2]):
    df = pd.read_excel(path, encoding = encoding, header=header)
    df.index = df[df.columns[2]]
    df = df.drop([df.columns[0], df.columns[1], df.columns[2]], axis=1)
    df.index.name = 'Freq [MHz]'
    return df

