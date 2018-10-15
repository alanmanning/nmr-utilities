import pylab as plt
plt.style.use('ggplot')
from IPython.core.debugger import Tracer; debug_here = Tracer()

def grid_plot(nrows,ncols,figsize,margins=(0.1,0.1,0.1,0.1),spaces=(0.04,0.04)):
    lmargin = margins[0]
    bmargin = margins[1]
    rmargin = margins[2]
    tmargin = margins[3]
    hspace = spaces[0]
    vspace = spaces[1]
    dx = (1-lmargin-rmargin-(ncols-1)*hspace)/ncols
    dy = (1-bmargin-tmargin-(nrows-1)*vspace)/nrows
    f = plt.figure(figsize=figsize)
    ax=[]
    for i in range(nrows):
        ax.append([])
        for j in range(ncols):
            pos = (lmargin+j*(hspace+dx),bmargin+(nrows-1-i)*(vspace+dy),dx,dy)
            ax[-1].append(f.add_axes(pos))

    return f,ax

def format_axes(f,ax_list,samex=True,samey=True,
        xlabel=None,ylabel=None,title=None,
        col_labels=None,row_labels=None):
    #figure out axes limits
    for i in list(range(len(ax_list))):
        for j in list(range(len(ax_list[0]))):
            if (i,j) == (0,0):
                max_xlim = list(ax_list[i][j].get_xlim())
                max_ylim = list(ax_list[i][j].get_ylim())
            else:
                xlim = ax_list[i][j].get_xlim()
                ylim = ax_list[i][j].get_ylim()                
                if xlim[0] < max_xlim[0]:
                    max_xlim[0] = xlim[0]
                if xlim[1] > max_xlim[1]:
                    max_xlim[1] = xlim[1]
                if ylim[0] < max_ylim[0]:
                    max_ylim[0] = ylim[0]
                if ylim[1] > max_ylim[1]:
                    max_ylim[1] = ylim[1]

    max_xlim[0]-=0.05*abs(max_xlim[0])
    max_xlim[1]+=0.05*abs(max_xlim[1])
    max_ylim[0]-=0.05*abs(max_ylim[0])
    max_ylim[1]+=0.05*abs(max_ylim[1])
    for i in list(range(len(ax_list))):
        for j in list(range(len(ax_list[0]))):
            #set the scales to be the same
            if samex:
                ax_list[i][j].set_xlim(max_xlim)
                if i != len(ax_list)-1:
                    ax_list[i][j].set_xticklabels([])
            if samey:
                ax_list[i][j].set_ylim(max_ylim)
                if j != 0:
                    ax_list[i][j].set_yticklabels([])

    if xlabel is not None:
        f.text(xlabel[0],xlabel[1],xlabel[2],
                {'ha':'center','va':'bottom','size':14})
    if ylabel is not None:
        f.text(ylabel[0],ylabel[1],ylabel[2],
                {'ha':'left','va':'center','size':14,'rotation':'vertical'})
    if title is not None:
        f.text(title[0],title[1],title[2],
                {'ha':'center','va':'bottom','size':12,'weight':'bold','color':'blue'})

    if col_labels is not None:
        for i in range(len(ax_list[0])):
            ax_list[0][i].set_title(col_labels[i],{'size':12,'weight':'bold','color':'black'})

    if row_labels is not None:
        for i in range(len(ax_list)):
            ax_list[i][0].set_ylabel(row_labels[i],{'size':12,'weight':'bold','color':'black'})

    return

def draw_legend(f,ax_list,ypos,ncol):
    for i in range(len(ax_list)):
        for j in range(len(ax_list[0])):
            handles,labels = ax_list[i][j].get_legend_handles_labels()
            if (handles,labels) != ([],[]):
                top_l = ax_list[0][0].get_position()
                top_r = ax_list[0][-1].get_position()
                pos = (top_l.xmin,ypos,top_r.xmax-top_l.xmin,(1-top_l.ymax)/2)
                leg_ax = f.add_axes(pos)
                legend = leg_ax.legend(handles,labels,ncol=ncol,frameon=True,loc='lower center',framealpha=1)
                legend.get_frame().set_facecolor('white')
                legend.get_frame().set_edgecolor('black')
                leg_ax.set_axis_off()
                return

def add_plot_text(ax,text,loc,fontsize):
    ax.text(loc[0], loc[1],text,
    horizontalalignment='left',
    verticalalignment='center',
    transform = ax.transAxes,
    fontsize = fontsize)
    return

def add_fig_title(fig,text):
    fig.suptitle(text,fontsize=12)
    return
