
from scipy import *
from pylab import *
from scipy.stats.kde import gaussian_kde
import scipy.stats.distributions
import scipy.stats

import pdb 

def get_coord(i, j, i_len):
    row = i + 1 
    col = j + 2 
    i_len += 2
    return row * i_len + col
    
def regression_t_test(x,y):
    beta_hat = dot(x,y)/dot(x,x) 
    alpha = mean(y) - beta_hat * mean(x)

    residuals = y - (alpha + beta_hat * x)
    sse = dot(residuals,residuals)
    x_hat = x - mean(x)

    t_score = beta_hat * sqrt(len(x)-2) / sqrt(sse / dot(x_hat,x_hat))

    return t_score

def scatter_grid(a, labels=None, mode="upper"):
    a = a.T
    fig = Figure()
    fig.set_facecolor('w')
    fig.set_edgecolor('w')

    for i in xrange(len(a)):
        for j in xrange(len(a)):
            if i >= j: continue

            if mode=="upper":
                sp = subplot(len(a)+2, len(a)+2, get_coord(i,j,len(a)), axisbg=(.9,.9,.9))
            else:
                sp = subplot(len(a)+2, len(a)+2, get_coord(j,i,len(a)), axisbg=(.9,.9,.9))
                
            sp.set_axis_off()
            if mode=="upper":
                x,y = j,i
            else:
                x,y = i,j
            scatter(a[x], a[y], zorder=1, color='k', marker='+')
            axhspan(0,100,0,1, facecolor=(.95,.95,.95), edgecolor=(.95,.95,.95), zorder=-100)
            sp.get_axes().set_aspect(1)

            m,b,r,p,se = scipy.stats.linregress(a[j],a[i])
            if p < 0.05:
                text(10,80,"p=%0.4f" %p,zorder=2, color=(.6,.6,.6))
                if mode=="upper":
                    plot((0,100), (b,b + m*100), color=(.8,.8,.8), zorder=.5)
                else:
                    plot((b,b + m*100), (0,100), color=(.8,.8,.8), zorder=.5)

            axis([0,100,0,100])
            
    for i in range(1,len(a)):
        if mode=="upper":
            sp = subplot(len(a)+2, len(a)+2, i+2)
        else:
            sp = subplot(len(a)+2, len(a)+2, (len(a)+1) * (len(a) +2) + i+1)
            
        sp.set_axis_off()

        if mode=="upper":
            i_h = i
            i_v = i-1
            bounds = 0,0.15
        else:
            i_h = i-1
            i_v = i
            bounds = 0.85, 1
        for p in a[i_h]:
            axvline(p, ymin=bounds[0], ymax=bounds[1], color='k')

        axis([0,100,0,100])

        k = gaussian_kde(a[i_h])
        points = arange(-100,300)
        estimate = k.evaluate(points)
        estimate /= sum(estimate)
        estimate = estimate[100:200]
        estimate = estimate * 1000+25
        if mode!='upper':
            estimate = 100 - estimate

        plot(points[100:200], estimate, color=(.6,.6,.6))
        axis([0,100,0,100])

        if mode=="upper":
            sp = subplot(len(a)+2, len(a)+2, get_coord(i-1, len(a), len(a)))
        else:
            sp = subplot(len(a)+2, len(a)+2, get_coord(i, -1, len(a)))
        sp.set_axis_off()

        for p in a[i_v]:
            axhline(p, xmin=bounds[0], xmax=bounds[1], color='k')
        axis([0,100,0,100])

        k = gaussian_kde(a[i_v])
        points = arange(-100,300)
        estimate = k.evaluate(points)
        estimate /= sum(estimate)
        estimate = estimate[100:200]
        estimate = estimate * 1000+25
        if mode!='upper':
            estimate = 100 - estimate

        plot(estimate, points[100:200], color=(.6,.6,.6))
        axis([0,100,0,100])

    if not labels is None:
        for i in xrange(len(a)):
            sp = subplot(len(a)+2,len(a)+2, (len(a)+2) * (i+1) + i+2)
            sp.set_axis_off()
            text(10,10,labels[i], rotation=45, color=(0.6,0.6,0.6), zorder=2)
            axis([0,100,0,100])

        
def full_scatter_grid(a,b,labels):
    scatter_grid(a, labels, mode='upper')
    scatter_grid(b, None, mode='lower')

def hash_and_density_plot(a, fig):
    x=0
    xwidth = 100./(len(a)+2)

    for i in xrange(0,len(a)):
        
        k = gaussian_kde(a[i])
        points = arange(-100,300)
        estimate = k.evaluate(points)
        estimate /= sum(estimate)
        estimate = estimate[100:200]
        estimate = estimate * 1000/(len(a)+1)+ x + xwidth/2

        plot(estimate+.5, points[100:200], color=(.6,.6,.6), linewidth=2)


        bounds = x,x+xwidth/2
        axvspan(x,x+xwidth/2,0,1, facecolor=(.95,.95,.95), edgecolor=(.95,.95,.95), zorder=-100)
        for p in a[i]:
            axhline(p, xmin=bounds[0]/100, xmax=bounds[1]/100, color='k', linewidth=2)

        
        #if not labels is None:
            #for i in xrange(len(a)):
                #sp = subplot(1, len(a),i+1)
                #xlabel(labels[i], color=(0.6,0.6,0.6), zorder=2)
    
        x+= 100./(len(a)+1)

    axis([0,100,0,100])

    fig.get_axes().set_axis_off()

def hash_and_density_plots(data_sets):
    
    category_count = data_sets[0].shape[1]
    
    for i in xrange(category_count):
        sp = subplot(1,category_count,i+1)
        
        attribute_data = [data_set[:,i] for data_set in data_sets]
        hash_and_density_plot(attribute_data, sp)

    subplots_adjust(wspace=0.)

def kspc_vs_wpm(k,w):
    scatter(k[:,0], w[:,0], marker=(5,2,0), color=(0.7,0.7,0.7), s=40)
    scatter(k[:,1], w[:,1], marker='^', color=(0.7,0.7,0.7), s=40)
    scatter(k[:,2], w[:,2], marker='s', color=(0.7,0.7,0.7), s=40)
    scatter(mean(k[:,0]), mean(w[:,0]), marker=(5,2,0), color='k', s=160)
    scatter(mean(k[:,1]), mean(w[:,1]), marker='^', color='k', s=160)
    scatter(mean(k[:,2]), mean(w[:,2]), marker='s', color='k', s=160)
    axis((1.,1.6,0,16))
