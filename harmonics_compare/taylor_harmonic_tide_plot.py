import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use('TkAgg')




class TaylorTide():
    def __init__(self,
                r_obs = None,
                rms_amp_max = 0.61,
                rms_amp_contours = [0.2, 0.4, 0.6],
                rms_err_contours=[0.2, 0.4, 0.6],
                cos_theta_lines = [0.3, 0.6, 0.9],
                theta_lines = [0, 15, 30, 45, 60],
                theta_lines_flag = True, # theta or cos(theta) construction lines
                err_contour_flag = True,
                ):

        self.fig, self.ax = self.plot_frame(r_obs, rms_amp_max, rms_amp_contours,
                                            rms_err_contours, cos_theta_lines, theta_lines,
                                            theta_lines_flag, err_contour_flag)

    # This custom formatter removes trailing zeros, e.g. "1.0" becomes "1", and
    def rms_fmt(self,x):
        s = f"{x:.2f}"
        if s.endswith("00"):  # ends with two zeros
            s = f"{x:.0f}"
        elif s.endswith("0"): # end with ONLY one zero
            s = f"{x:.1f}"
        return rf"{s}"

    def plot_frame(self, r_obs, rms_amp_max, rms_amp_contours, rms_err_contours,
                   cos_theta_lines, theta_lines,
                   theta_lines_flag, err_contour_flag):

        #fig = plt.figure(figsize=(2, 2)) # to make a thumbnail schematic
        fig = plt.figure()

        theta = np.arange(0, np.pi + np.pi/100, np.pi/100)

        # define meshgrid for contour plots
        x = np.arange(0,rms_amp_max, rms_amp_max/100)
        y = np.arange(0,rms_amp_max, rms_amp_max/100)
        X,Y = np.meshgrid(x,y)

        # setting the axis limits in [left, bottom, width, height]
        rect = [0.1, 0.1, 0.8, 0.8]

        # the cartesian axis:
        ax = fig.add_axes(rect, frameon=False)
        #ax =fig.add_subplot(111)

        # RMS amplitude arc contours
        if rms_amp_contours != []:
            Camp = ax.contour( X,Y,np.sqrt(X**2 + Y**2), levels=rms_amp_contours, colors='grey', linestyles='dotted')
            ax.clabel(Camp, Camp.levels, inline=True, fmt=self.rms_fmt, fontsize=10)

        # Obs point, arc through obs, RMS error arcs
        if r_obs is not None:
            ax.scatter(r_obs, 0, s=30, color='blue', clip_on=False)  # obs point

            if err_contour_flag is True:
                ax.plot(r_obs * np.cos(theta), r_obs * np.sin(theta), '-', color='blue')  # arc through obs

                # RMS error arc from obs as origin
                mask = X**2 + Y**2 > rms_amp_max**2
                C = np.ma.masked_where(mask, np.sqrt((X-r_obs)**2 + Y**2))
                Cerr = ax.contour(X, Y, C, levels=rms_err_contours, colors='grey',
                                  linestyles='dashed')
                ax.clabel(Cerr, Cerr.levels, inline=True, fmt=self.rms_fmt, fontsize=10)

                # Add text contour label. THIS WILL BE A PROBLEM IF MORE THAN 3 LEVELS ARE DEFINED
                fmt = {}
                strs = ['', '', 'RMS error']
                for l, s in zip(rms_err_contours, strs):
                    fmt[l] = s
                ax.clabel(Cerr, Cerr.levels, inline=True, fmt=fmt, fontsize=10)


        # Bounding lines - black
        ax.plot( rms_amp_max*np.cos(theta), rms_amp_max*np.sin(theta), '-', color='black')
        # Bounding x=0 and y=0 lines - black
        ax.plot( [0,rms_amp_max], [0,0], '-', color='black')
        ax.plot( [0,0], [0,rms_amp_max], '-', color='black')


        # Cos theta / correlation lines
        r = np.arange(0, rms_amp_max, rms_amp_max/50)
        if theta_lines_flag == False:
            for cos_theta in cos_theta_lines:
                ax.plot( r*cos_theta, r*np.sqrt(1 - cos_theta**2), ':', color='grey')
                ax.text( rms_amp_max*cos_theta, rms_amp_max*np.sqrt(1 - cos_theta**2), str(cos_theta), color='k' )
            ax.text( rms_amp_max*1/np.sqrt(2), rms_amp_max*1/np.sqrt(2), "Correlation", rotation=-45, color='k')
        else:
            for ang in theta_lines:
                ang_rad = ang*np.pi/180.
                ax.plot( r*np.cos(ang_rad), r*np.sin(ang_rad), ':', color='grey')
                ax.text( rms_amp_max*np.cos(ang_rad), rms_amp_max*np.sin(ang_rad), str(ang)+"$^o$", color='k' )
            ax.text( rms_amp_max*1/np.sqrt(2), rms_amp_max*1/np.sqrt(2), "phase error", rotation=-45, color='k')


        # axis limits and axis labels
        ax.set_xlim([0,rms_amp_max])
        ax.set_ylim([0,rms_amp_max])
        ax.set_xlabel('RMS amplitude (m)')
        ax.set_ylabel('RMS amplitude (m)')

        ax.set_aspect(1)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if(0):
            # the polar axis:
            ax_polar = fig.add_axes([0.1, 0.0, 0.8, 0.8], polar=True, frameon=False)
            #ax_polar = fig.add_axes(ax.get_position(original=True), polar=True, frameon=False)
            ## [left, bottom, width, height]
            #ax_polar.set_position([0.1, 0.0, 0.8, 0.8])

            # the polar plot
            ax_polar.plot(r, r, color='r', linewidth=3)
            ax_polar.set_rmax(2.0)
            ax_polar.grid(True)
            ax_polar.set_thetamin(0)
            ax_polar.set_thetamax(54)
            ax_polar.set_xlabel('RMS amplitude (m)')
            ax_polar.set_ylabel('RMS amplitude (m)')


        return fig, ax

if __name__ == '__main__':
    # Add data
    R = ['1.00', '1.00', '0.81', '0.84', '1.00', '1.00', '0.90', '0.91']
    rms_amp= ['0.51', '0.51', '0.63', '0.51', '0.41', '0.40', '0.45', '0.38']
    rms_err= ['0.00', '0.02', '0.37', '0.29', '0.00', '0.01', '0.20', '0.17']
    label = ['obs:s', 'fes:s', 'gs1p1:s', 'gs1p2:s', 'obs:d', 'fes:d', 'gs1p1:d', 'gs1p2:d']

    # put the data into arrays
    rms_amp = np.array([float(x) for x in rms_amp])
    rms_err = np.array([float(x) for x in rms_err])
    R = np.array([float(x) for x in R])

    # Create TaylorTide plot template
    tt = TaylorTide(
        r_obs=rms_amp[0],
        rms_amp_max=0.7,
        rms_amp_contours=[0.2, 0.4, 0.6],
        rms_err_contours=[0.2, 0.4, 0.6],
        cos_theta_lines=[0.3, 0.6, 0.9],
        )
    # Add data to axes
    tt.ax.scatter( rms_amp[1:4] * R[1:4], rms_amp[1:4] * np.sqrt(1 - R[1:4]**2) , s=10, c='r')

    # Create TaylorTide plot template
    tt = TaylorTide(
        r_obs = rms_amp[4],
        rms_amp_max=0.61,
        rms_amp_contours=[0.2, 0.4, 0.6],
        rms_err_contours=[0.2, 0.4, 0.6],
        cos_theta_lines=[0.3, 0.6, 0.9],
        )
    # Add data to axes
    tt.ax.scatter( rms_amp[5::] * R[5::], rms_amp[5::] * np.sqrt(1 - R[5::]**2) , s=10, c='r')

    tt = TaylorTide(
        r_obs=0.4,
        rms_amp_max=0.6,
        rms_amp_contours=[],
        rms_err_contours=[0.2, 0.4],
        theta_lines=[0, 30, 60],
        theta_lines_flag=True,  # theta or cos(theta) construction lines
        err_contour_flag=True,
        )
    # Add data to axes
    #tt.ax.scatter( rms_amp[5::] * R[5::], rms_amp[5::] * np.sqrt(1 - R[5::]**2) , s=10, c='r')
    plt.show()