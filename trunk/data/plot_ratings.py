
from scipy import *
from pylab import *
from scipy.stats.kde import gaussian_kde
import scipy.stats.distributions
import scipy.stats

import pdb 

def get_coord(i, j, i_len):
    row = i + 1 
    col = j + 1 
    i_len += 1
    return row * i_len + col
    
def regression_t_test(x,y):
    beta_hat = dot(x,y)/dot(x,x) 
    alpha = mean(y) - beta_hat * mean(x)

    residuals = y - (alpha + beta_hat * x)
    sse = dot(residuals,residuals)
    x_hat = x - mean(x)

    t_score = beta_hat * sqrt(len(x)-2) / sqrt(sse / dot(x_hat,x_hat))

    return t_score


def scatter_grid(a, labels=None):
    a = a.T
    fig = Figure()
    fig.set_facecolor('w')
    fig.set_edgecolor('w')

    for i in xrange(len(a)):
        for j in xrange(len(a)):
            if i >= j: continue
            
            sp = subplot(len(a)+1, len(a)+1, get_coord(i,j,len(a)), axisbg=(.9,.9,.9))
            sp.set_axis_off()
            scatter(a[j], a[i], zorder=1, color='k', marker='+')
            axhspan(0,100,0,1, facecolor=(.95,.95,.95), edgecolor=(.95,.95,.95), zorder=-100)
            sp.get_axes().set_aspect(1)

            m,b,r,p,se = scipy.stats.linregress(a[j],a[i])
            if p < 0.01:
                text(10,80,"p=%0.4f" %p,zorder=2, color=(.6,.6,.6))
                plot((0,100), (b,b + m*100), color=(.6,.6,.6))

            axis([0,100,0,100])
            
    for i in range(1,len(a)):
        sp = subplot(len(a)+1, len(a)+1, i+1)
        sp.set_axis_off()

        for p in a[i]:
            axvline(p, ymin=0, ymax=0.15, color='k')
        axis([0,100,0,100])

        k = gaussian_kde(a[i])
        points = arange(-100,300)
        estimate = k.evaluate(points)
        estimate /= sum(estimate)
        estimate = estimate[100:200]
        plot(points[100:200], estimate*1000+25, color='k')
        axis([0,100,0,100])

        sp = subplot(len(a)+1, len(a)+1, get_coord(i-1, len(a), len(a)))
        sp.set_axis_off()

        for p in a[i-1]:
            axhline(p, xmin=0, xmax=0.15, color='k')
        axis([0,100,0,100])

        k = gaussian_kde(a[i-1])
        points = arange(-100,300)
        estimate = k.evaluate(points)
        estimate /= sum(estimate)
        estimate = estimate[100:200]
        plot(estimate*1000+25, points[100:200], color='k')
        axis([0,100,0,100])

    if not labels is None:
        for i in xrange(len(a)):
            sp = subplot(len(a)+1,len(a)+1, (len(a)+1) * (i+1) + i+1)
            sp.set_axis_off()
            text(10,10,labels[i], rotation=45, color=(0.6,0.6,0.6), zorder=2)
            axis([0,100,0,100])

        

