import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import dbp_general as dbg
#import scipy.stats as stats

def harmonic_polar_plot_auto_loc(des_locs, mod_locs, obs_locs):
    
    return

def harmonic_polar_plot(a_mod, g_mod, a_obs, g_obs, 
                        const = ['M2','S2','K1','K2','N2','P1','O1','Q1'], 
                        title = '', cmplx_sum = False): 
    '''
    # For comparison of individual and total harmonics at specific individual 
    # locations. Will plot up all harmonics in (r, theta). Modelled harmonics 
    # will be plotted with a solid line and observed harmonics with a dashed 
    # line.
    #
 IN # a_mod      :: vector of model amplitudes for each constituent.
 IN # g_mod      :: vector of model phases for each constituent.
 IN # a_obs      :: vector of observed amplitudes.
 IN # g_obs      :: vector of observed phases.
 IN # const      :: constituent order of data rows (used for subtitles)
 IN # title      :: title for the plot.
 IN # cmplx_sum  :: whether or not to plot a vector sum of all harmonic
    #
OUT # None
    '''
    
    plt.figure(figsize = (6,6))
    if cmplx_sum == True:
        const.append('Sum')
    
    #Plot model lines
    multi_polar_plot(a_mod, g_mod, style = '-', cmplx_sum = cmplx_sum)
    #Legend
    leg = plt.legend(const,loc=9,bbox_to_anchor=(0.5,0), ncol=5, framealpha=1)
    leg.get_frame().set_linewidth(0)
    #Plot observed lines
    multi_polar_plot(a_obs, g_obs, style = '--', cmplx_sum = cmplx_sum)
        
    #Central dot for aesthetic reasons
    plt.polar([0,0], color=[0.5,0.5,0.5], marker = 'o', markersize = 3.5)
    
    #Title
    plt.title(title, loc='right')
    
    return

def multi_polar_plot(r, theta, cmap=cm.get_cmap('tab10'), style = '-', 
                     title = '', degrees = True, cmplx_sum = False):
    '''
    # For putting multiple polar plots on one figure, with option of adding a 
    # summed polar plot (done in complex space). 
    #
 IN # r         :: Vector of radii
 IN # theta     :: vector of angles
 IN # degrees   :: True if the input is in degrees (is converted to radians)
 IN # cmplx_sum :: If true then an additional line is plotted which is the 
    #              sum of the inputs in complex space.
    #
OUT # fig, ax   :: Matplotlib figure/axis objects for the plot
    '''
    
    #Check lengths of input arrays.
    n_points = len(r)
    n_colors = cmap.N
    
    if degrees:
        theta = np.radians(theta)

    #Loop over each dataset and point and plot a line      
    for pp in range(0,n_points):
        color_tmp = cmap(pp%n_colors)
        ax = plt.polar([theta[pp], theta[pp]], [0, r[pp]], color = color_tmp, linewidth = 1.25,
                        marker = '.', markersize=3, linestyle = style)
        
    if cmplx_sum:
        z1, z2 = dbg.polar2cart(r,theta)
        z1_sum = np.sum(z1); z2_sum = np.sum(z2)
        r_sum, theta_sum = dbg.cart2polar(z1_sum, z2_sum)
        plt.polar([theta_sum, theta_sum], [0, r_sum], color = 'k', linewidth = 1.25,
                  marker = '.', markersize = 3, linestyle = style)
    
    return plt.gcf(), ax

def qqplot(X, Y, c = [1,0.7,0.4], quantiles = 'auto', yex = True, 
           subplots=(1,1), title='', subtitles=[''], xl = '', yl = ''):
    '''
    # Quantile-quantile plot. Will plot multiple subplots if  specified.
    # In this case the rows of the input data are used for each individual 
    # subplot. By default, all data is plotted on a single subplot. Also 
    # removes NaNs from the data. Data is be plotted from left to right, top 
    # to bottom.
    #
 IN # X,Y       :: [Necessary] Input data
 IN # c         :: Fill color of markers [default orangey]
 IN # quantiles :: List of quantile boundaries. [default 0.25 quantiles]
 IN # yex       :: If true then plot y = x line. [default True]
 IN # subtitles :: List of titles for each individual subplot (if using). Will 
    #           :: loop if len(list) != number of plots [default none]
 IN # xl, yl    :: x and y label for entire plot (not subplots) [default none]
 IN # subplots  :: (rows, cols) for array of subplots. [default (1,1)]
    #
OUT # fig, ax   :: Matplotlib figure/axis objects for plot.
    '''
    
    if type(quantiles) == str:
        quantiles = np.arange(0.05,1,0.025)
    
    #General dimensions and shapes of input variables
    input_shape = np.shape(X)
    n_dims = len(input_shape)
    nr, nc = subplots
    
    # Check the form of the input array: 2D or 1D? If 1D (of the form (n,))
    # expand the dimensions so that it works with later code.
    if n_dims == 2:
        n_plots, n_dp = input_shape # Number of plots and number of datapoints
    elif n_dims == 1:
        n_plots = 1
        X = np.expand_dims(X,1); X = np.transpose(X)
        Y = np.expand_dims(Y,1); Y = np.transpose(Y)
    else:
        raise Exception('Input arrays of wrong dimension')
   
    # Define figure properties, overall title and axis labels
    f, ax = plt.subplots(nr, nc, figsize=(nc*3+1.5,nr*3+1))
    f.suptitle(title)
    f.text(0.5, 0.04, xl, ha='center')
    f.text(0.04, 0.5, yl, va='center', rotation='vertical')
    
    # Loop over each subplot, from left to right, top to bottom
    for plot_ii in range(0,n_plots):
        
        # Get relevant line from input data and remove nans
        X_tmp = X[plot_ii,:]
        Y_tmp = Y[plot_ii,:]
        idx = np.isfinite(X_tmp) & np.isfinite(Y_tmp)
        X_tmp = X_tmp[idx]
        Y_tmp = Y_tmp[idx]
        
        #Calculate quantiles for this plot.
        qq_x = np.quantile(X_tmp, quantiles)
        qq_y = np.quantile(Y_tmp, quantiles)
        
        # Define current subplot axis
        if nr == 1 and nc == 1:
            ax_tmp = ax
        elif nr == 1 or nc == 1:
            ax_tmp = ax[int(plot_ii)]
        else:
            # This calculation is to ensure the left -> right, top -> bottom
            # scheme
            ax_tmp = ax[int(plot_ii/nc), plot_ii-nc*int(plot_ii/nc)]
    
        # Define axis limits as the maximum(xmax, ymax) + 10% (square plots)
        # These are also used for plotting lines.
        Xmin = np.nanmin(qq_x); Xmax = np.nanmax(qq_x)
        Ymin = np.nanmin(qq_y); Ymax = np.nanmax(qq_y)
        axmax = np.max([Xmax,Ymax]); axmin = np.min([Xmin,Ymin])
        axrange10 = 0.1*(axmax - axmin)
        
        # If y=x line is to be plotted then plot it
        if yex:
            lineX = [axmin,axmax]
            fityx = np.poly1d([1,0])
            ax_tmp.plot(lineX, fityx(lineX), c=[0.5,0.5,0.5], linewidth=1.5)

        #Plot quantiles
        ax_tmp.scatter(qq_x, qq_y, marker='o', edgecolor='k', color = c, 
                       alpha=0.8)

        # Set axis limits, set square axes and add subplot 
        # title.
        ax_tmp.set_xlim(axmin-axrange10, axmax+axrange10)
        ax_tmp.set_ylim(axmin-axrange10, axmax+axrange10)
        ax_tmp.set_aspect('equal', adjustable='box')
        # Set individual subtitles
        n_subs = len(subtitles)
        ax_tmp.title.set_text(subtitles[plot_ii%n_subs])
                  
    return f, ax

def scatter_with_fit(X, Y, C=None, dofit = True, yex = True, subplots=(1,1),
                     title='', subtitles=[''], xl = '', yl = '', qualcmap = 0,
                     qualbounds = [], cmap='cmo.deep'):
    '''
    # Scatter plot with a linear fit to the data. Can plot multiple subplots if 
    # specified. In this case the rows of the input data are used for each 
    # individual subplot. Will also remove NaNs from the data. Data will be 
    # plotted from left to right, top to bottom.
    #
 IN # X,Y       :: [Necessary] Input data
 IN # C         :: Vector of colour data or a single colour for all points.
    #              [default none]
 IN # dofit     :: If true then plot linear fit. [default true]
 IN # yex       :: If true then plot y = x line. [default true]
 IN # subtitles :: List of titles for each individual subplot (if using). Will 
    #           :: loop if len(list) != number of plots [default none]
 IN # xl, yl    :: x and y label for entire plot (not subplots) [default none]
 IN # subplots  :: (rows, cols) for array of subplots. [default (1,1)]
 IN # qualcmap  :: If true, use a qualitative colormap defined by qualbounds.
    #           :: [default false]
    #
OUT # fig, ax   :: Matplotlib figure/axis objects for plot.
    '''
    
    #General dimensions and shapes of input variables
    input_shape = np.shape(X)
    n_dims = len(input_shape)
    nr, nc = subplots
    
    # Check the form of the input array: 2D or 1D? If 1D (of the form (n,))
    # expand the dimensions so that it works with later code.
    if n_dims == 2:
        n_plots, n_dp = input_shape # Number of plots and number of datapoints
    elif n_dims == 1:
        n_plots = 1
        X = np.expand_dims(X,1); X = np.transpose(X)
        Y = np.expand_dims(Y,1); Y = np.transpose(Y)
    else:
        raise Exception('Input arrays of wrong dimension')
        
    # Load cmap if necessary
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
   
    # Define figure properties, overall title and axis labels
    f, ax = plt.subplots(nr, nc, figsize=(nc*3+1.5,nr*3+1))
    f.suptitle(title)
    f.text(0.5, 0.04, xl, ha='center')
    f.text(0.04, 0.5, yl, va='center', rotation='vertical')
    
    # Loop over each subplot, from left to right, top to bottom
    for plot_ii in range(0,n_plots):
        
        # Get relevant line from input data and remove nans
        X_tmp = X[plot_ii,:]
        Y_tmp = Y[plot_ii,:]
        idx = np.isfinite(X_tmp) & np.isfinite(Y_tmp)
        X_tmp = X_tmp[idx]
        Y_tmp = Y_tmp[idx]
        
        # Define current subplot axis
        if nr == 1 and nc == 1:
            ax_tmp = ax
        elif nr == 1 or nc == 1:
            ax_tmp = ax[int(plot_ii)]
        else:
            # This calculation is to ensure the left -> right, top -> bottom
            # scheme
            ax_tmp = ax[int(plot_ii/nc), plot_ii-nc*int(plot_ii/nc)]
    
        # Define axis limits as the maximum(xmax, ymax) + 10% (square plots)
        # These are also used for plotting lines.
        Xmin = np.nanmin(X_tmp); Xmax = np.nanmax(X_tmp)
        Ymin = np.nanmin(Y_tmp); Ymax = np.nanmax(Y_tmp)
        axmax = np.max([Xmax,Ymax]); axmin = np.min([Xmin,Ymin])
        axrange10 = 0.1*(axmax - axmin)
        
        # If y=x line is to be plotted then plot it
        if yex:
            lineX = [axmin,axmax]
            fityx = np.poly1d([1,0])
            ax_tmp.plot(lineX, fityx(lineX),c=[0.5,0.5,0.5],linewidth=1.5)
            
        # If fit is to be plotted then fit it and plot it
        if dofit:
            lineX = [Xmin,Xmax]
             #Calculate data fit and cast to poly1d object
            fit_tmp = np.polyfit(X_tmp, Y_tmp, 1)
            fit = np.poly1d( fit_tmp )
            ax_tmp.plot(lineX, fit(lineX),c=[1,128/255,0],linewidth=1.5)
            rs = dbg.r_squared_lin(X_tmp, Y_tmp, fit)
            
        # COLORS: Define scatter colormap if necessary (if C has been defined)
        # If no color vector defined then set to black
        if C is None:
            im = ax_tmp.scatter(X_tmp, Y_tmp ,marker='o', s=10, c='k',
                                edgecolor='none')
        else:
            # Remove indices where there were NaNs from the C array
            C_tmp = C[idx]
            if qualcmap:
                # If qualitative colormap is to be used, define the boundaries.
                # cmap example: colors.ListedColormap([[0, 204/255, 0], 
                #               [1, 128/255, 0], [0,102/255,204/255]])                    
                boundaries = qualbounds #e.g. [0, 100, 300, 10000]
                norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
                im = ax_tmp.scatter(X_tmp, Y_tmp ,marker='o', s=10, c=C_tmp, 
                                    norm=norm, cmap=cmap, edgecolor='none')
            else:
                #Continuous colormap
                im = ax_tmp.scatter(X_tmp, Y_tmp ,marker='o', s=10, c=C_tmp,
                                    edgecolor='none', cmap=cmap)
            # Plot Colorbar on a new subplot to the right 
            f.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.25, hspace=0.1)
            cb_ax = f.add_axes([0.83, 0.15, 0.01, 0.7])
            f.colorbar(im, cax=cb_ax, extend='max')

        # Set axis limits, set square axes and add subplot 
        # title.
        ax_tmp.set_xlim(axmin-axrange10, axmax+axrange10)
        ax_tmp.set_ylim(axmin-axrange10, axmax+axrange10)
        ax_tmp.set_aspect('equal', adjustable='box')
        # Set individual subtitles
        n_subs = len(subtitles)
        ax_tmp.title.set_text(subtitles[plot_ii%n_subs])
        
        #fit text
        if dofit:
            ax_tmp.text(0.4,0.125,'{} {:03.2f} {} {:03.2f}'.format('y =',
                        fit_tmp[0],'x +',fit_tmp[1]),
                        transform=ax_tmp.transAxes)
            ax_tmp.text(0.4,0.05,'{} {:03.2f} '.format('$R^2$ =',
                        rs),transform=ax_tmp.transAxes)
                
    return f, ax

