#%%```python
from common_variables import *

# Define paths
fd = './Figures/manuscripts/'
os.makedirs(fd, exist_ok=True)
indir = './data/'

def create_figure_1B():
    """Create and save Figure 1B"""
    # Load data
    fn = './data/Figure_1B.data.pkl'
    with open(fn, 'rb') as f:
        data = pickle.load(f)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(4, 4))
    
    # Create outer pie chart
    outer_wedge, texts = ax.pie(data['odata'],
                               wedgeprops=dict(width=0.3, antialiased=True),
                               colors=data['colors'],
                               startangle=90,
                               shadow=False,
                               textprops=data['textprops'])
    
    # Create inner pie chart
    inner_wedge, texts = ax.pie(data['idata'],
                               wedgeprops=dict(width=0.3, antialiased=True),
                               labels=data['ilbs'],
                               colors=data['colors'],
                               startangle=90,
                               shadow=False,
                               textprops=data['textprops'],
                               labeldistance=0.7,
                               radius=0.68)
    
    # Add titles
    ax.set_title('normal', fontsize=12, y=1.)
    ax.text(0, 0, 'tumor', fontsize=12, va='center', ha='center')
    
    # Configure annotation settings
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props,
              zorder=0,
              va="center")
    
    # Add annotations for outer wedges
    for i, p in enumerate(outer_wedge):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(data['olbs'][i],
                   xy=(x, y),
                   xytext=(1.1 * np.sign(x), 1.2 * y),
                   horizontalalignment=horizontalalignment,
                   **kw)
    return fig

def create_figure_1C():
    """Create and save Figure 1C"""
    # Load data
    fn = './data/Figure_1C.data.tsv.gz'
    df = pd.read_csv(fn, sep='\t')
    df['cohort'] = ['Tumor' if x == 'Cancer' else 'Normal' for x in df['cohort']]
    tissues = np.unique(df['primary_diagnosis'])

    fig, axs = plt.subplots(figsize=(9, 4), ncols=3, nrows=2)
    fig.subplots_adjust(wspace=0.4)
    axs = axs.flatten()
    for i, tissue in enumerate(tissues):
        D = df[df['primary_diagnosis'] == tissue]
        D = D.sort_values(['feature', 'cohort'], ascending=[False, True])
        significance = []
        grouped_vals = D.groupby(['feature', 'cohort'])['log2Enrich'].apply(list)
        pvals = []
        for feature in D['feature'].unique():
            stat, pval = sp.stats.ranksums(*(grouped_vals[feature].values))
            pvals.append(pval)
        fdr = sm.stats.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        annot = []
        for __ in pvals:
            if __ < 0.0001:
                annot.append('****')
            elif __ < 0.001:
                annot.append('***')
            elif __ < 0.01:
                annot.append('**')
            elif __ < 0.05:
                annot.append('*')
            else:
                annot.append('')
        for k, text in enumerate(annot):
            axs[i].text(2.75, k+0.25, text, fontsize=8, ha='center', va='center');
        
        # Fix boxplot to show boxes by using the x and y parameters correctly
        sns.boxplot(x='log2Enrich', y='feature', hue='cohort', data=D, 
                    hue_order=['Normal', 'Tumor'], 
                    palette=cmaps['cohort'], 
                    fliersize=0, 
                    ax=axs[i], 
                    linewidth=1,
                    showbox=True,
                    width=0.7);
        
        axs[i].set_xlim(-2.5, 3.1);
        axs[i].set_ylabel('');
        axs[i].set_xlabel('');
        axs[i].tick_params(labelsize=8)
        axs[i].text(s=tissue, x=0.05, y=0.95, va='top', ha='left', fontsize=10, transform=axs[i].transAxes, weight='bold')
        if i != 4:
            axs[i].legend().remove()
        else:
            axs[i].legend(frameon=False, bbox_to_anchor=(1, 1), fontsize=8)
    axs[-1].set_visible(False);
    axs[-2].set_xlabel(r'Enrichment (log$_2$)');
    axs[-3].set_xlabel(r'Enrichment (log$_2$)');
    return fig

def create_figure_2A():
    """Create and save Figure 2A"""
    fn = os.path.join(indir, 'Figure_2.tissue_drg.per_tissue.data.pkl')
    with open(fn, 'rb') as f:
        tissue_drg_yh = pickle.load(f)

    groups = sorted(list(tissue_drg_yh.keys()))
    fig, axs = plt.subplots(figsize=(6, 6), ncols=2, nrows=3, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.03, hspace=0.2)
    axs = axs.flatten()

    for i, group in enumerate(groups):
        D = tissue_drg_yh[group]
        D['kind'] = 'not significant'
        D = D.sort_values('FDR', ascending=False)
        D.loc[(D['FDR'] < 0.05)&(D['FC'] < 0), 'kind'] = 'down-regulated'
        D.loc[(D['FDR'] < 0.05)&(D['FC'] > 0), 'kind'] = 'up-regulated'
        lbs, cnt = (np.unique(D['kind'], return_counts=True))
        dict_obj  = {x:y for x,y in zip(lbs, cnt)}
        D['color'] = D['kind'].map(colors)

        axs[i].scatter(D['logCPM'], D['logFC'], c=D['color'].values, s=0.5);

        axs[i].set_title('{:}'.format(group), y=0.97);
        axs[i].text(0.9, 0.8, '{:}'.format(dict_obj['up-regulated']), transform=axs[i].transAxes, ha='right', color=colors['up-regulated'], weight='bold');
        axs[i].text(0.95, 0.05, '{:}'.format(dict_obj['down-regulated']), transform=axs[i].transAxes, ha='right', color=colors['down-regulated'], weight='bold');
    axs[4].set_xlabel(r'log$_2$(CPM)', fontsize=12);
    axs[5].set_xlabel(r'log$_2$(CPM)', fontsize=12);
    axs[2].set_ylabel(r'log$_2$(tumor/normal FC)', fontsize=12);
    return fig

def create_figure_2B():
    """Create and save Figure 2B"""
    fn = os.path.join(indir, 'Figure_2.tissue_drg.per_tissue.data.pkl')
    with open(fn, 'rb') as f:
        tissue_drg_yh = pickle.load(f)

    tissue_drg_pl = pd.concat([tissue_drg_yh[x]['logFC'].rename(x) for x in tissue_drg_yh], axis=1)
    cols = ['breast', 'colon', 'lung', 'ovary', 'pancreas', 'all']
    tissue_drg_pl = tissue_drg_pl[cols]

    fig, axs = plt.subplots(figsize=(10, 8.5), ncols=6, nrows=6)
    fig.subplots_adjust(wspace=0, hspace=0);
    corr_axes = []

    xlims_x = pd.DataFrame(index=range(len(cols)), columns=range(len(cols)))
    xlims_y = pd.DataFrame(index=range(len(cols)), columns=range(len(cols)))
    ylims_x = pd.DataFrame(index=range(len(cols)), columns=range(len(cols)))
    ylims_y = pd.DataFrame(index=range(len(cols)), columns=range(len(cols)))
    lin_regs = pd.DataFrame(index=range(len(cols)), columns=range(len(cols)))

    for i in range(len(cols)):
        for j in range(1+i, len(cols)):
            ax = axs[i, j]
            ax.set_visible(False);

    for i in range(1, len(cols)):
        for j in range(0, i):
            c1 = cols[i]
            c2 = cols[j]
            ax = axs[i, j]
            x = tissue_drg_pl[c1]
            y = tissue_drg_pl[c2]
            lm = sm.OLS(endog=y, exog=sm.add_constant(x)).fit()
            _b, _m = lm.params
            rsq = round(lm.rsquared_adj, 3)
            density2d(x=x, y=y, ax=ax, s=2)
            if j != 0:
                ax.set_yticks([])
            else:
                ax.set_yticks([-2, -1, 0, 1, 2])
            if i != len(cols)-1:
                ax.set_xticks([])
            else:
                ax.set_xticks([-2, -1, 0, 1, 2])
            ax.plot(x, _b + x*_m, c='k', lw=1, alpha=1)
            text = 'y = {:.2f} + {:.2f} x\n'.format(_b, _m)+r'R$^2$: '+'{:.3f}'.format(rsq)
            ax.text(0.05, 0.8, text, transform=ax.transAxes, fontsize=9, ha='left', linespacing=0.8)
            corr_axes.append(ax)
            _const = 1
            xlims_x.loc[i, j] = min(x)*_const
            xlims_y.loc[i, j] = max(x)*_const
            ylims_x.loc[i, j] = min(y)*_const
            ylims_y.loc[i, j] = max(y)*_const

    xlims = []
    for i, j in zip(xlims_x.min(0), xlims_y.max(0)):
        xlims.append((i, j))
    ylims = []
    for i, j in zip(ylims_x.min(1), ylims_y.max(1)):
        ylims.append((i, j))

    for i in range(1, len(cols)):
        ylim = ylims[i]
        axs[i, 0].set_ylabel(cols[i])
        for j in range(0, i):
            xlim = xlims[j]
            ax = axs[i, j]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if i == len(cols) - 1:
                axs[i, j].set_xlabel(cols[j])

    for i, key in enumerate(cols):
        ax = axs[i, i]
        ax.hist(tissue_drg_pl[key], bins=50, density=True);
        ax.text(0.9, 0.8, key, transform=ax.transAxes, ha='right');
        if i != len(cols) - 1:
            ax.set_xticks([])
        ax.set_yticks([])
        xlim = xlims[i]
        ax.set_xlim(xlim)

    return fig
    
def create_figure_2E():
    """Create and save Figure 2E"""
    fn = os.path.join(indir, 'Figure_2E.pathway_NES.tsv.gz')
    D = pd.read_csv(fn, sep='\t', index_col=0)
        
    gfg = nhm(data=D, figsize=(14, 7), cmapCenter='coolwarm', linewidths=0, xrot=90, tick_size=10, cmaps={'center':0})
    gfg.hcluster(row_cluster=True, col_cluster=False)
    fig, axs = gfg.run(ax_gap=10, center_args={'mid':0, 'max':4, 'min':-4})
    axs[0][1].tick_params(axis='x', labelsize=12)
    yticklbs = axs[0][1].get_yticklabels()
    axs[0][3].set_title('NES')
    new_yticks = []

    return gfg.fig

def create_figure_3A():
    tumor_tissue_early_late_drg = pickle.load(open(os.path.join(indir, 'Figure_3A.tumor_tissue_early_late_drg.pkl'), 'rb'))
    cfdna_early_drg = pickle.load(open(os.path.join(indir, 'Figure_3A.cfdna_early_drg.pkl'), 'rb'))
    in_genes = np.intersect1d(tumor_tissue_early_late_drg['breast'].index, cfdna_early_drg['breast'].index)

    groups = sorted(list(tumor_tissue_early_late_drg.keys()))
    # groups = sorted(list(tumor_tissue_early_late_drg.keys()))[1:]
    tumor_tissue_drg = {}
    tumor_tissue_models = {}
    # fig, axs = plt.subplots(figsize=(9, 6), ncols=3, nrows=2)
    fig, axs = plt.subplots(figsize=(7, 4.5), ncols=3, nrows=2)
    fig.subplots_adjust(wspace=0.15, hspace=0.25)
    axs = axs.flatten()

    for i, group in enumerate(groups):
        D = tumor_tissue_early_late_drg[group].loc[in_genes]
        tumor_tissue_drg[group] = D
        density2d(data=D, x='early.logFC', y='late.logFC', ax=axs[i], bins=500, xlabel='', ylabel='', s=1);
        # add_identity(axs[i], ls='--', color='red');

        model = sm.OLS(D['late.logFC'], sm.add_constant(D['early.logFC'])).fit()
        intercept, slope = model.params
        intercept = round(intercept, 2)
        slope = round(slope, 2)
        tumor_tissue_models[group] = model

        density2d(data=D, x='early.logFC', y='late.logFC', ax=axs[i], bins=500, xlabel='', ylabel='', s=3);
        r2 = round(model.rsquared, 2)
        slope_pval = round(model.pvalues[1], 4)
        f_pval = round(model.f_pvalue, 2)

        t_test = model.t_test('early.logFC = 1')
        p_slope_equals_one = round(t_test.pvalue.item(), 4)

        axs[i].plot(D['early.logFC'], D['early.logFC']*slope + intercept, color='k', lw=1, alpha=.75)
        # axs[i].text(0.05, 0.95, r'y = {:} + {:} x'.format(intercept, slope) + '\n' + r'R$^2$: ' + '{:}\n'.format(r2) + r'$p < 0.05$', ha='left', va='top', transform=axs[i].transAxes, fontsize=8)
        axs[i].text(0.05, 0.95, r'y = {:} + {:} x'.format(intercept, slope) + '\n' + r'R$^2$: ' + '{:}\n'.format(r2), ha='left', va='top', transform=axs[i].transAxes, fontsize=10)
    
        axs[i].set_title('{:}'.format(group), y=1);
    axs[1].text(0.5, 1.15, 'tumor tissues', ha='center', transform=axs[1].transAxes, fontsize=12);
    axs[0].text(-0.4, 0.5, r'late', fontsize=12, ha='left', va='center', transform=axs[0].transAxes, rotation=90);
    axs[3].text(-0.4, 0.5, r'late', fontsize=12, ha='left', va='center', transform=axs[3].transAxes, rotation=90);
    axs[-2].set_xlabel(r'early', fontsize=12, ha='center');
    
    return fig

def create_figure_3C():
    D = pd.read_csv(os.path.join(indir, 'Figure_3C.tissues.tumor_fraction_by_stage.tsv.gz'), sep='\t', index_col=0)
    significance = []
    grouped_vals = D.groupby(['tissue', 'early_late_stage'])['tumor_fraction'].apply(list)
    pvals = []
    for tissue in D['tissue'].unique():
        stat, pval = sp.stats.ranksums(*(grouped_vals[tissue].values))
        pvals.append(pval)
    fdr = sm.stats.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
    annot = []
    for __ in fdr:
        annot.append('*' if __ < 0.05 else 'n.s.')   


    fig, ax = plt.subplots(figsize=(8, 2))
    colors = np.concatenate([['C{:}'.format(x), 'C{:}'.format(x)] for x in range(0, 5)])
    sns.boxplot(data=D, x='tissue', y='tumor_fraction', ax=ax, hue='early_late_stage', dodge=True, fliersize=0, color='darkred')
    np.random.seed(42)
    sns.stripplot(data=D, x='tissue', y='tumor_fraction', ax=ax, hue='early_late_stage', dodge=True, jitter=True, color='black', alpha=0.7, legend=False);
    ax.set_xticklabels([x.get_text().split(' - ')[-1] for x in ax.get_xticklabels()], rotation=0)
    ax.legend(bbox_to_anchor=(1, 1), frameon=False, title='stage');
    ax.set_xlabel('cancer type', fontsize=10);
    ax.set_ylabel('tumor fraction');
    ax.tick_params(labelsize=10, axis='x', pad=12);
    ax.set_ylim(-0.05, 1.05);
    ax.set_title('tumor tissues', fontsize=10);
    for i, text in enumerate(annot):
        if text == '*':
            ax.text(i, 0.9, text, fontsize=12, ha='center', va='center', transform=ax.transData);
        else:
            ax.text(i, 0.9, text, fontsize=10, ha='center', va='center', transform=ax.transData);
        ax.text(i-0.3, -0.14, 'early', fontsize=10)
        ax.text(i+0.11, -0.14, 'late', fontsize=10)

    return fig

def create_figure_4A():
    """Create and save Figure 4A"""
    drg = pickle.load(open(os.path.join(indir, 'Figure_4A.cfDNA_drg.pkl'), 'rb'))
    groups = sorted(list(drg.keys()))
    fig, axs = plt.subplots(figsize=(8, 4), ncols=3, nrows=2, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.03, hspace=0.2)
    axs = axs.flatten()

    for i, group in enumerate(groups):
        D = drg[group].sort_values('sig')
        lbs, cnt = (np.unique(D['kind'], return_counts=True))
        dict_obj  = {x:y for x,y in zip(lbs, cnt)}
        color_d = D.set_index('kind')['color'].to_dict()

        axs[i].scatter(D['logCPM'], D['logFC'], c=D['color'].values, s=0.5);

        axs[i].set_title('{:}'.format(group), y=0.97);
        axs[i].text(0.9, 0.8, '{:}'.format(dict_obj['up-regulated']), transform=axs[i].transAxes, ha='right', color=color_d['up-regulated'], weight='bold');
        axs[i].text(0.95, 0.05, '{:}'.format(dict_obj['down-regulated']), transform=axs[i].transAxes, ha='right', color=color_d['down-regulated'], weight='bold');
    axs[4].set_xlabel(r'log$_2$(CPM)', fontsize=12);
    axs[0].set_ylabel(r'log$_2$(cancer/control FC)', fontsize=12, y=0);
    
    return fig

def create_figure_4B():
    """Create and save Figure 4B"""
    D = pd.read_csv(os.path.join(indir, 'Figure_4B.cfDNA.tumor_fraction_by_stage.tsv.gz'), sep='\t', index_col=0)
    grouped_vals = D.groupby(['tissue', 'early_late_stage'])['tumor_fraction'].apply(list)
    pvals = []
    for tissue in D['tissue'].unique():
        stat, pval = sp.stats.ranksums(*(grouped_vals[tissue].values))
        pvals.append(pval)
    fdr = sm.stats.multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
    annot = []
    for __ in fdr:
        annot.append('*' if __ < 0.05 else 'n.s.')
    fig, ax = plt.subplots(figsize=(8, 2))
    colors = np.concatenate([['C{:}'.format(x), 'C{:}'.format(x)] for x in range(0, 5)])
    sns.boxplot(data=D, x='tissue', y='tumor_fraction', ax=ax, hue='early_late_stage', dodge=True, fliersize=0, color='darkred')
    np.random.seed(42)
    sns.stripplot(data=D, x='tissue', y='tumor_fraction', ax=ax, hue='early_late_stage', dodge=True, jitter=True, color='black', alpha=0.7, legend=False);
    ax.set_xticklabels([x.get_text().split(' - ')[-1] for x in ax.get_xticklabels()], rotation=0)
    ax.legend(bbox_to_anchor=(1, 1), frameon=False, title='stage');
    ax.set_xlabel('cancer type', fontsize=10);
    ax.set_ylabel('tumor fraction');
    ax.tick_params(labelsize=10, axis='x', pad=12);
    ax.set_ylim(-0.05, 1.05);
    ax.set_title('cancer cfDNA', fontsize=10);
    for i, text in enumerate(annot):
        if text == '*':
            ax.text(i, 0.9, text, fontsize=12, ha='center', va='center', transform=ax.transData);
        else:
            ax.text(i, 0.9, text, fontsize=10, ha='center', va='center', transform=ax.transData);
        ax.text(i-0.3, -0.14, 'early', fontsize=10)
        ax.text(i+0.11, -0.14, 'late', fontsize=10)

    return fig

def create_figure_4C():
    tumor_tissue_early_late_drg = pickle.load(open(os.path.join(indir, 'Figure_3A.tumor_tissue_early_late_drg.pkl'), 'rb'))
    cfdna_early_drg = pickle.load(open(os.path.join(indir, 'Figure_3A.cfdna_early_drg.pkl'), 'rb'))
    cfdna_late_drg = pickle.load(open(os.path.join(indir, 'Figure_3A.cfdna_late_drg.pkl'), 'rb'))
    in_genes = np.intersect1d(tumor_tissue_early_late_drg['breast'].index, cfdna_early_drg['breast'].index)

    groups = sorted(list(cfdna_early_drg.keys()))
    cfdna_drg = {}
    cfdna_models = {}
    fig, axs = plt.subplots(figsize=(7, 4.5), ncols=3, nrows=2)
    fig.subplots_adjust(wspace=0.2, hspace=0.3)
    axs = axs.flatten()

    for i, group in enumerate(groups):
        D_early = cfdna_early_drg[group]
        D_late = cfdna_late_drg[group]
        D = pd.DataFrame({'early.logFC':D_early.loc[in_genes, 'logFC'],
                        'late.logFC':D_late.loc[in_genes, 'logFC'],
                        'early.FDR':D_early.loc[in_genes, 'FDR'],
                        'late.FDR':D_late.loc[in_genes, 'FDR']},
                        index=in_genes)
        cfdna_drg[group] = D
        # add_identity(axs[i], ls='--', color='red');
        model = sm.OLS(D['late.logFC'], sm.add_constant(D['early.logFC'])).fit()
        cfdna_models[group] = model
        intercept, slope = model.params
        intercept = round(intercept, 2)
        slope = round(slope, 2)

        density2d(data=D, x='early.logFC', y='late.logFC', ax=axs[i], bins=500, xlabel='', ylabel='', s=3);
        r2 = round(model.rsquared, 2)
        slope_pval = round(model.pvalues[1], 4)
        f_pval = round(model.f_pvalue, 2)

        t_test = model.t_test('early.logFC = 1')
        p_slope_equals_one = round(t_test.pvalue.item(), 4)

        axs[i].plot(D['early.logFC'], D['early.logFC']*slope + intercept, color='k', lw=1, alpha=.75)
        # axs[i].text(0.05, 0.95, r'y = {:} + {:} x'.format(intercept, slope) + '\n' + r'R$^2$: ' + '{:}\n'.format(r2) + r'$p < 0.05$', ha='left', va='top', transform=axs[i].transAxes, fontsize=8)
        axs[i].text(0.05, 0.95, r'y = {:} + {:} x'.format(intercept, slope) + '\n' + r'R$^2$: ' + '{:}\n'.format(r2), ha='left', va='top', transform=axs[i].transAxes, fontsize=10)
        axs[i].set_title('{:}'.format(group), y=0.97);
    axs[1].text(0.5, 1.15, 'cfDNA', ha='center', transform=axs[1].transAxes, fontsize=12);
    axs[0].text(-0.4, 0.5, r'late', fontsize=12, ha='left', va='center', transform=axs[0].transAxes, rotation=90);
    axs[3].text(-0.4, 0.5, r'late', fontsize=12, ha='left', va='center', transform=axs[3].transAxes, rotation=90);
    axs[-2].set_xlabel(r'early', fontsize=12, ha='center');

    return fig

def create_figure_4D():
    """Create and save Figure 4D"""
    res = pd.read_csv(os.path.join(indir, 'Figure_4D.cfDNA_concordance_by_tissue.tsv.gz'), sep='\t', index_col=0)
    chi2_p = []
    for key in res.index:
        x1 = res.loc[key, 'early concordant up'] + res.loc[key, 'early concordant down']
        xp = (res.loc[key, 'early total significant'] - x1)
        y1 = res.loc[key, 'late concordant up'] + res.loc[key, 'late concordant down']
        yp = (res.loc[key, 'late total significant'] - y1)
        cont_tab = np.array([[x1, xp], [y1, yp]]).astype(int)
        chi2_p.append(sp.stats.chi2_contingency(cont_tab)[1])
    res['chi2_pval'] = chi2_p
    threshold = 0.05
    D = res[['early % concordant', 'late % concordant']]
    D.columns = ['early (I+II)', 'late (III+IV)']

    fig, ax = plt.subplots(figsize=(5, 3))
    D.plot(kind='bar', ax=ax, rot=0, legend=True, linewidth=1, color=['snow', 'darkred'], edgecolor='black', width=0.75)
    legend = ax.get_legend()
    legend.set_title('stage')

    # Get the maximum y-value
    ymax = D.max(1)

    # Iterate over each xtick
    for i, xtick in enumerate(ax.get_xticks()):
        barplot_annotate_brackets(ax, xtick, ymax[i], '*', dx=.2)

    ax.set_ylim(0, 30);
    ax.set_ylabel('% concordant DhmR genes')
    ax.set_title('tumor tissue and cancer cfDNA')
    ax.legend(frameon=False);
    return fig
    
def create_figure_4E():
    """Create and save Figure 4E"""
    annotate_pathways = [
        'RAF_UP.V1_DN',
        'BRCA1_DN.V1_UP',
        'LEF1_UP.V1_DN',
        'KRAS.600_UP.V1_DN',
        'NFE2L2.V2',
        'KRAS.BREAST_UP.V1_UP',
        'KRAS.50_UP.V1_UP',
        'KRAS.300_UP.V1_UP',
        'KRAS.LUNG.BREAST_UP.V1_UP',
        'KRAS.600_UP.V1_UP',
        'KRAS.600.LUNG.BREAST_UP.V1_UP']
    D = pd.read_csv(os.path.join(indir, 'Figure_4E.cfDNA_tissue_concordant_GSEA.tsv.gz'), sep='\t', index_col=0)
    dfc = pd.DataFrame(D.columns, columns=['tissue'])
    gfg = nhm(data=D,
            figsize=(18, 15),
            cmapCenter='coolwarm',
            linewidths=0,
            xrot=90,
            tick_size=8,
            showyticks=True,
            cmaps={
                'center': 0,
            })
    gfg.hcluster(row_cluster=True, col_cluster=False)
    fig, axs = gfg.run(ax_gap=10, center_args={'mid': 0, 'max': 4, 'min': -4})
    axs[0][1].tick_params(axis='x', labelsize=12)
    axs[0][3].set_title('NES')
    yticklbs = axs[0][1].get_yticklabels()
    for tklb in yticklbs:
        txt = tklb.get_text()
        if txt in annotate_pathways:
            tklb.set_weight('bold')
    return fig

def create_figure_5A():
    """Create and save Figure 5A"""
    D = pd.read_csv(os.path.join(indir, 'Figure_5A.UMAP.tsv.gz'), sep='\t', index_col=0)
    s = 10
    fig, axs = plt.subplots(figsize=(7, 3), ncols=2)

    for ax in axs:
        for spine in ax.spines:
            ax.spines[spine].set_visible(False)
        ax.set_yticks([])
        ax.set_xticks([]);

    axs[0].scatter(D['umap_x'], D['umap_y'], c=D['tissue_normal_color'], s=s, alpha=0.8);
    axs[0].set_title('normal')
    axs[0].set_xlabel('UMAP 1')
    axs[0].set_ylabel('UMAP 2')
    # for i, tissue in enumerate(D['tissue_tumor_color'].unique()):
    for i, tissue in enumerate(np.unique(D['tissue_tumor'])):
        if tissue != 'None':
            axs[1].scatter(D[D['tissue_tumor'] == tissue]['umap_x'],
                            D[D['tissue_tumor'] == tissue]['umap_y'],
                            c=D[D['tissue_tumor'] == tissue]['tissue_tumor_color'],
                            label = tissue,
                            s=s, alpha=0.8);
        else:
            axs[1].scatter(D[D['tissue_tumor'] == tissue]['umap_x'],
                            D[D['tissue_tumor'] == tissue]['umap_y'],
                            c=D[D['tissue_tumor'] == tissue]['tissue_tumor_color'],
                            s=s, alpha=0.8);
    axs[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=8);
    axs[1].set_title('tumor');
    axs[1].set_xlabel('UMAP 1')
    axs[1].set_ylabel('UMAP 2');
    return fig
    
def create_figure_5B():
    """Create and save Figure 5B"""
    D = pickle.load(open(os.path.join(indir, 'Figure_5B.model_score.pkl'), 'rb'))
    df = D['df']
    df.index = [x.split('.')[0] for x in df.index]
    dfr = D['dfr']
    dfc = D['dfc']
    colormaps = {'cohort': {'Cancer': 'red', 'Normal': 'blue'},
                'tissue': cmaps['primary_diagnosis']}

    fig, axs = nhm(data=df,
                dfr=dfr,
                dfc=dfc,
                wspace=0.005,
                srot=90,
                figsize=(25, 5),
                cmaps={
                    'tissue': {'breast':'C0', 'colon':'C1', 'lung':'C2', 'ovary':'C3', 'pancreas':'C4'},
                    'cohort': {'Cancer':'red', 'Normal':'blue'}
                }).run(center_args={'cbar_title':'prediction score'})
    return fig
    
def create_figure_5C():
    """Create and save Figure 5C"""
    roc_curves_to_plot = pd.read_csv(os.path.join(indir, 'Figure_5C.roc_curve.cfdna_toto.tsv.gz'), sep='\t', index_col=0)
    palette = roc_curves_to_plot.set_index('label')['color'].to_dict()

    fig, ax = plt.subplots(figsize=(5, 3))
    sns.lineplot(data=roc_curves_to_plot, x='1 - specificity', y='sensitivity', hue='label', palette=palette, ax=ax);
    ax.legend(frameon=True, fontsize=8);
    ax.set_title('Tissue of Tumor Origin Prediction in cfDNA samples', fontsize=10);
    add_identity(ax, xlim=(0, 1), ylim=(0, 1));
    return fig
    
def create_figure_5D():
    """Create and save Figure 5D"""
    df_top_two = pd.read_csv(os.path.join(indir, 'Figure_5D.confusion_matrix_top_two.cfdna_toto.tsv.gz'), sep='\t', index_col=0)
    sensitivity = np.diag(df_top_two.values) / df_top_two.sum(0)
    xticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(sensitivity.index, sensitivity.values)]
    specificity = calculate_specificity(df_top_two)
    yticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(specificity.index, specificity.values)][::-1]

    fig, ax = plt.subplots(figsize=(5, 4.5))
    gfg = sns.heatmap(df_top_two.iloc[::-1], annot=True, cmap='Blues', fmt='d', cbar=False, ax=ax, annot_kws={'fontsize': 12}, vmin=0, vmax=50);
    ax.set_xlabel('True Class (sensitivity)', fontsize=12, labelpad=8);
    ax.set_ylabel('Predicted Class (specificity)', fontsize=12, labelpad=8);
    ax.set_xticklabels(xticks, ha='center', fontsize=10);
    ax.set_yticklabels(yticks, ha='center', fontsize=10, x=-0.03);
    return fig

def create_figure_5E():
    """Create and save Figure 5E"""
    roc_curves_to_plot = pd.read_csv(os.path.join(indir, 'Figure_5E.roc_curve.cfdna_toto_on_tumor_tissues.tsv.gz'), sep='\t', index_col=0)
    palette = roc_curves_to_plot.set_index('label')['color'].to_dict()

    fig, ax = plt.subplots(figsize=(5, 3))
    sns.lineplot(data=roc_curves_to_plot, x='1 - specificity', y='sensitivity', hue='label', palette=palette, ax=ax);
    ax.legend(frameon=True, fontsize=8);
    ax.set_title('Tissue of Tumor Origin Prediction in tumor tissue samples', fontsize=10);
    add_identity(ax, xlim=(0, 1), ylim=(0, 1));
    return fig

def create_figure_5F():
    """Create and save Figure 5F"""
    df_top_two = pd.read_csv(os.path.join(indir, 'Figure_5F.confusion_matrix_top_two.cfdna_toto_on_tumor_tissues.tsv.gz'), sep='\t', index_col=0)
    sensitivity = np.diag(df_top_two.values) / df_top_two.sum(0)
    xticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(sensitivity.index, sensitivity.values)]
    specificity = calculate_specificity(df_top_two)
    yticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(specificity.index, specificity.values)][::-1]

    fig, ax = plt.subplots(figsize=(5, 4.5))
    gfg = sns.heatmap(df_top_two.iloc[::-1], annot=True, cmap='Blues', fmt='d', cbar=False, ax=ax, annot_kws={'fontsize': 12}, vmin=0, vmax=50);
    ax.set_xlabel('True Class (sensitivity)', fontsize=12, labelpad=8);
    ax.set_ylabel('Predicted Class (specificity)', fontsize=12, labelpad=8);
    ax.set_xticklabels(xticks, ha='center', fontsize=10);
    ax.set_yticklabels(yticks, ha='center', fontsize=10, x=-0.03);
    return fig


if __name__ == "__main__":
    # Create Figure 1B
    fig_1B = create_figure_1B() 
    savefig(fig_1B, fd+'Figure1B')
    
    # Create Figure 1C
    fig_1C = create_figure_1C()
    savefig(fig_1C, fd+'Figure1C')

    # Create Figure 2A
    fig_2A = create_figure_2A() 
    savefig(fig_2A, fd+'Figure2A')

    # Create Figure 2B
    fig_2B = create_figure_2B()
    savefig(fig_2B, fd+'Figure2B')

    # Create Figure 2D
    fig_2E = create_figure_2E()
    savefig(fig_2E, fd+'Figure2E')
   
    # Create Figure 3A
    fig_3A = create_figure_3A()
    savefig(fig_3A, fd+'Figure3A')
   
    # Create Figure 3C
    fig_3C = create_figure_3C()
    savefig(fig_3C, fd+'Figure3C')
   
    # Create Figure 4A
    fig_4A = create_figure_4A()
    savefig(fig_4A, fd+'Figure4A')

    # Create Figure 4B
    fig_4B = create_figure_4B()
    savefig(fig_4B, fd+'Figure4B')

    # Create Figure 4C
    fig_4C = create_figure_4C()
    savefig(fig_4C, fd+'Figure4C')

    # Create Figure 4D
    fig_4D = create_figure_4D()
    savefig(fig_4D, fd+'Figure4D')

    # Create Figure 4E
    fig_4E = create_figure_4E()
    savefig(fig_4E, fd+'Figure4E')

    # Create Figure 5A
    fig_5A = create_figure_5A()
    savefig(fig_5A, fd+'Figure5A')

    # Create Figure 5B
    fig_5B = create_figure_5B()
    savefig(fig_5B, fd+'Figure5B')

    # Create Figure 5C
    fig_5C = create_figure_5C()
    savefig(fig_5C, fd+'Figure5C')

    # Create Figure 5D
    fig_5D = create_figure_5D()
    savefig(fig_5D, fd+'Figure5D')

    # Create Figure 5E
    fig_5E = create_figure_5E()
    savefig(fig_5E, fd+'Figure5E')

    # Create Figure 5F
    fig_5F = create_figure_5F()
    savefig(fig_5F, fd+'Figure5F')