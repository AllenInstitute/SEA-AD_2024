import os
import copy
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import interpolate as sp_interpolate
from scipy import stats as sp_stats
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_fscore_support as scores
from statsmodels.gam.api import GLMGam, BSplines
from patsy.builtins import *
import warnings

warnings.filterwarnings("ignore")
pwd = os.getcwd()

def smooth(x, y, xgrid):
    samples = np.random.choice(len(x), int(len(x) * 0.8), replace=True)
    y_s = y[samples]
    x_s = x[samples]
    y_sm = sm_lowess(y_s, x_s, return_sorted = False)
    # regularly sample it onto the grid
    y_grid = sp_interpolate.interp1d(x_s, y_sm, fill_value="extrapolate")(xgrid)
    return y_grid

def movingaverage_nan(d, window_size):
    new_d = np.zeros_like(d)
    d= np.array(d, dtype=np.float64)
    for i in np.arange(d.shape[0]):
        interval = [np.max([0, i - int(np.round(window_size / 2))]), np.min([d.shape[0] - 1, i + int(np.round(window_size / 2))])]
        value = np.nanmean(d[interval[0]:interval[1]])
        new_d[i] = value
        
    return new_d

def fit_gam(df, outcome, spline_var, formula, n_gauss=5, deg_free=3, alpha=1e-3):
    np.seterr(invalid='ignore')
    try:
        x_spline = df[spline_var]
        bs = BSplines(x_spline, df=n_gauss + 1, degree=deg_free)
        gam_bs = GLMGam.from_formula(
            formula = "Q('" + outcome + "') ~ " + formula,
            data=df,
            smoother=bs,
            alpha=alpha,
        )
        results = gam_bs.fit()
        
        good_return = []
        for i in range(n_gauss):
            good_return.extend([results.params[spline_var + "_s" + str(i)]])
            
        for i in range(n_gauss):
            good_return.extend([results.pvalues[spline_var + "_s" + str(i)]])
            
        good_return.extend([results.pseudo_rsquared()])
        
        return good_return
    
    except:
        failed_return = np.zeros((1, n_gauss)).tolist()[0]
        failed_return.extend(np.zeros((1, n_gauss)).tolist()[0])
        failed_return.extend([0])
        return failed_return
    

def generate_results(snames, fname):
    out = pd.DataFrame()
    for s_ in sel_names:
        model_outputs = fit_gam(df=qn, outcome=s_, spline_var=spline_var, formula=formula)
        out = pd.concat([out, pd.DataFrame(model_outputs).T], axis = 0)

    out.index = sel_names
    out.columns = ['spline0', 'spline1', 'spline2', 'spline3', 'spline4', 'pval0', 'pval1', 'pval2', 'pval3', 'pval4', 'r2']

    out
    plt.figure()
    sns.heatmap(-out[out.columns[out.columns.str.contains('spline')]], cmap="Blues")
    plt.figure()
    sns.heatmap(-out[out.columns[out.columns.str.contains('pval')]], cmap="Blues")
    plt.savefig(fname)
    
    return out
    

def preprocess_smooth(x_slot, y_slot, pand):
    x = pand[x_slot].to_numpy()
    y = pand[y_slot].to_numpy()
    xgrid = np.linspace(x.min(),x.max())
    K = 1000
    smooths = np.stack([smooth(x, y, xgrid) for k in range(K)]).T   
    mean = np.nanmean(smooths, axis=1)
    stderr = scipy.stats.sem(smooths, axis=1)
    stderr = np.nanstd(smooths, axis=1, ddof=0)
    
    return xgrid, mean, stderr

def lmplots(
        data,
        feature_name,
        y_scale="relative",
        to_plot="Cell type",
        to_plot_filter=None,
        figsize=None,
        dpi=100,
        celltype_filter=None,
        cmap="tab10",
        level_order=None):

    # y scale transformations
    if y_scale == "relative":
        sample_sums = np.sum(data.X, axis=1, keepdims=True)
        X = data.X/sample_sums
        value_name = "Proportion"
    elif y_scale == "log_relative":
        sample_sums = np.sum(data.X, axis=1, keepdims=True)
        X = np.log10((data.X + 1)/(sample_sums + data.shape[1]))
        value_name = "log(Proportion)"
    elif y_scale == "log":
        X = np.log10(data.X + 1)
        value_name = "log(count)"
    elif y_scale == "count":
        X = data.X
        value_name = "count"
    else:
        raise ValueError("Invalid y_scale transformation")
        
    count_df = pd.DataFrame(X, columns=data.var.index, index=data.obs.index).\
        merge(data.obs[feature_name], left_index=True, right_index=True)

    plot_df = pd.melt(count_df, id_vars=feature_name, var_name="Cell type", value_name=value_name)
    if to_plot != "Cell type":
        plot_df = plot_df.\
            merge(data.var[to_plot], left_on="Cell type", right_index=True)

    
    if to_plot_filter is not None:
        plot_df = plot_df[plot_df[to_plot].isin(to_plot_filter)]
        
    if celltype_filter is not None:
        plot_df = plot_df[plot_df["Cell type"].isin(celltype_filter)]
    
    plt.rcParams["figure.figsize"] = figsize
    plot_df[to_plot] = plot_df[to_plot].cat.remove_unused_categories()
    for i,j in enumerate(plot_df[to_plot].cat.categories):
        data = plot_df[plot_df[to_plot] == j]
        x = data[feature_name].to_numpy()
        y = data[value_name].to_numpy()
        xgrid = np.linspace(x.min(),x.max())
        K = 1000
        smooths = np.stack([smooth(x, y, xgrid) for k in range(K)]).T   
        mean = np.nanmean(smooths, axis=1)
        stderr = sp_stats.sem(smooths, axis=1)
        stderr = np.nanstd(smooths, axis=1, ddof=0)
        # plot it
        plt.fill_between(xgrid, mean-2*stderr, mean+2*stderr, color=cmap[j], alpha=0.06)
        plt.plot(xgrid, mean, color=cmap[j], label=j)
        
    plt.tight_layout()
    
    return plt.gca()

def delta_plot(
    adata,
    genes,
    groupby,
    plotby,
    donor,
    across,
    title,
    groupby_subset=None,
    highlight=[],
    colormap=None,
    figsize=None,
    nrow=1,
    ncol=None,
    min_cells=500,
    save=False,
    normalize_to_start=False,
    legend=False,
    xlim=(0,1),
    ylim=None,
):
    vars_to_get = copy.copy(genes)
    vars_to_get.extend([groupby, plotby, across, donor])
    vars_to_get = np.unique(vars_to_get).tolist()
    tmp = sc.get.obs_df(adata, vars_to_get)
    tmp[across] = tmp[across].fillna(0)
    cell_counts = tmp[plotby].value_counts(sort=False) > min_cells
    cell_counts = cell_counts[cell_counts == True].index.to_list()
    tmp = tmp.loc[
        [j in cell_counts for j in tmp[plotby]],
        :
    ]
    tmp[title] = tmp.loc[:, genes].mean(axis=1)
    if title not in genes:
        tmp = tmp.drop(genes, axis=1)
        
    plot_categories = adata.obs[groupby].cat.categories
    if groupby_subset != None:
        plot_categories = groupby_subset
        
    if ncol == None:
        ncol = np.int32(len(plot_categories) / nrow)

    if figsize == None and legend == False:
        figsize = (3 * ncol, 5 * nrow)
    if figsize == None and legend == True:
        figsize = (4 * ncol, 5 * nrow)
    
    plt.rcParams["figure.figsize"] = figsize
    
    tmp[groupby] = tmp[groupby].astype("str")
    tmp[plotby] = tmp[plotby].astype("str")
    tmp[groupby + plotby] = tmp[groupby] + "|" + tmp[plotby]
    tmp[groupby + plotby] = tmp[groupby + plotby].astype("category")
    
    tmp = tmp.groupby([groupby, plotby, donor, groupby + plotby]).mean().dropna().reset_index()

    for c,d in enumerate(plot_categories):
        for l in tmp[groupby + plotby].cat.categories:
            if tmp.loc[tmp[groupby + plotby] == l, groupby].iloc[0] == d:
                line_kws = {"lw": 1, "alpha": 1}
                fill_kws = {"alpha": 0.1}
                ci = None
                scatter = False
                if highlight == [] and colormap is not None:
                    supertype_color = colormap[l.replace(d + "|", "")]
                else:
                    if l.replace(d + "|", "") in highlight:
                        supertype_color = "dodgerblue"
                    else:
                        supertype_color = "goldenrod"
            else:
                supertype_color = "lightgrey"
                line_kws = {"lw": 0.5, "alpha": 0.5}
                fill_kws = {"alpha": 0.001}
                ci = None
                scatter = False
                
            ax = plt.subplot(nrow, ncol, c+1);
            
            data = tmp.loc[(tmp[groupby + plotby] == l), :]
            x = data[across].to_numpy()
            y = data[title].to_numpy()
            xgrid = np.linspace(x.min(),x.max())
            K = 1000
            smooths = np.stack([smooth(x, y, xgrid) for k in range(K)]).T
            mean = np.nanmean(smooths, axis=1)
            stderr = np.nanstd(smooths, axis=1)
            # plot it
            plt.plot(
                xgrid,
                mean,
                color=supertype_color,
                label=l,
                linewidth=line_kws["lw"],
                alpha=line_kws["alpha"],
            )
            plt.fill_between(
                xgrid,
                mean-stderr,
                mean+stderr,
                color=supertype_color,
                alpha=fill_kws["alpha"]
            )
            
            plt.xlim(xlim);
            if ylim != None:
                plt.ylim(ylim);
            ax.set(title=d, xlabel=None, ylabel=None);
            
        if highlight != [] and colormap == None and legend == True and c == ncol - 1:
            leg = plt.legend(loc="upper left", labels=["No highlight", "Highlight"], bbox_to_anchor=(1.05, 1), title=plotby);
            for g,h in enumerate(leg.legendHandles):
                h.set_color(["goldenrod", "dodgerblue"][g])
                h.set_linewidth(5)
    
    if highlight == [] and colormap is not None and legend == True:
        leg = plt.legend(loc="center left", labels=adata.obs[plotby].cat.categories, bbox_to_anchor=(1.05, 0.5), ncol=2, title=plotby);
        for c,d in enumerate(leg.legendHandles):
            d.set_color(list(colormap.values())[c])
            d.set_linewidth(5)
    plt.suptitle(title);

    if normalize_to_start == True:
        plot_to_return = plt.gcf();
        for i in range(len(plot_to_return.axes)):
            sub_value = 0
            done = False
            for j in range(len(plot_to_return.axes[i].lines)):
                x,y = plot_to_return.axes[i].lines[j].get_data()
                sub_value = copy.copy(y[0])
                y = y - sub_value
                plot_to_return.axes[i].lines[j].set_data(x, y)
                
                path = plot_to_return.axes[i].collections[j].get_paths()[0]
                path._vertices = path._vertices - [0, sub_value]
                plot_to_return.axes[i].collections[j]._paths = [path]

        for i in range(len(plot_to_return.axes)):
            plot_to_return.axes[i].relim()
    else:
        plot_to_return = plt.gcf();
    plt.tight_layout(w_pad=1.2, h_pad=1.1)
    if save != False:
        warnings.filterwarnings("ignore")
        plt.savefig(save.replace("{title}", title), bbox_inches="tight");
    return plot_to_return

def plot_confusion(df, figsize=(20,20)):
    label_list = df['true'].cat.categories.to_list()
    C = confusion_matrix(df['true'],
                        df['pred'],
                        labels=label_list,
                        normalize='true')
    C_df = pd.DataFrame(C, index=label_list, columns=label_list)
    C_df.index.name = 'True'
    C_df.columns.name = 'Pred'
    plt.subplots(1,1,figsize=figsize)
    ax = sns.heatmap(C_df,
                vmin=0,vmax=1,
                xticklabels=True, 
                yticklabels=True, 
                square=True, 
                linewidths=0.1,
                cmap="Greys",
                cbar_kws={'shrink': 0.3,"orientation": "vertical"})

    return C_df

def get_scores(df):
    p, r, f1, s = scores(df['true'], df['pred'],
                         labels=df['true'].cat.categories.to_list(),
                         average=None,
                         sample_weight=None)

    result = pd.DataFrame({'label': df['true'].cat.categories,
                           'precision': p,
                           'recall': r,
                           'f1': f1,
                           'support': s})
    return result