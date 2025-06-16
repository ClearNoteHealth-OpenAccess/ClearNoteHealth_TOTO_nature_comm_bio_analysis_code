#%%```python
from common_variables import *

# Define paths
fd = './Figures/manuscripts/'
os.makedirs(fd, exist_ok=True)
indir = './data/'

def create_sfigure_1A():
    """Create and save Supplemental Figure 1A"""
    # Load data
    fn = './data/SFigure_1A.nheatmap.pkl'
    with open(fn, 'rb') as f:
        data = pickle.load(f)

    dist = data['dist']
    dfr = data['dfr']
    cmaps = data['cmaps']
    gfg = nhm(dist,
            dfr=dfr,
            figsize=(10.5, 5),
            xrot=90,
            linewidths=0,
            cmapCenter='gist_heat_r',
            srot=90,
            cmaps={'has_cancer':cmaps['has_cancer'],
                    'primary_diagnosis':{'breast':'C0', 'colon':'C1', 'lung':'C2', 'ovary':'C3', 'pancreas':'C4'},
                    'tissue_type':cmaps['tissue_type']},
            )
    gfg.hcluster(metric='correlation', method='average')
    fig, ax = gfg.run(center_args={'max':1, 'cbar_title':'correlation distance'})

    return fig

def create_sfigure_2A():
    """Create and save Supplemental Figure 2A"""
    drg_tissue_feature_type = pd.read_csv(os.path.join(indir, 'SFigure_2A.barplot.tsv.gz'), sep='\t', index_col=0)

    fig, ax = plt.subplots(figsize=(8, 2))
    drg_tissue_feature_type.plot(kind='bar', color=cmaps['primary_diagnosis'], ax=ax);
    ax.set_yscale('log');
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter())
    ax.set_title('percent differential features (FDR < 0.05)');
    ax.legend(bbox_to_anchor=(1, 1));
    ax.tick_params(axis='x', rotation=0);
    return fig
    
def create_sfigure_3A():
    """Create and save Supplemental Figure 3A"""
    fn = os.path.join(indir, 'SFigure_3A.boxplot.extracted_cfDNA_conc.tsv.gz')
    D = pd.read_csv(fn, sep='\t', index_col=0)
    corder = ['control', 'breast', 'colon', 'lung', 'ovary', 'pancreas'] 

    fig, ax = plt.subplots(figsize=(7.5, 2.5))
    sns.boxplot(D, x='tissue', y='dna_conc_ng_ul', hue='early_late_stage', palette='Reds', ax=ax,
                order=corder,
                hue_order=['control', 'early (I+II)', 'late (III+IV)']);
    ax.set_ylabel('concentration (ng/µL)')
    ax.set_title('extracted cfDNA')
    ax.set_xlabel('cancer type')
    ymax = D.groupby(['tissue'])['dna_conc_ng_ul'].max().loc[corder]
    for i in range(1, 6):
        barplot_annotate_brackets(ax, i+0.1, ymax[i], '*', bar_pad=10)
    ax.legend(bbox_to_anchor=(1, 1), frameon=False, title='stage');
    ax.set_ylim(0.18, 200);
    ax.set_yscale('log')
    return fig

def create_sfigure_3B():
    """Create and save Supplemental Figure 3B"""
    with open(os.path.join(indir, 'SFigure_3B.venn_diagram.pkl'), 'rb') as f:
        data = pickle.load(f)
    set1_up = data['set1_up']
    set2_up = data['set2_up']
    set1_down = data['set1_down']
    set2_down = data['set2_down']

    fig, axs = plt.subplots(figsize=(6, 3), ncols=2)
    fig.subplots_adjust(wspace=0)

    venn2([set1_up, set2_up], ('Tissue', 'cfDNA'), ax=axs[0]);
    axs[0].set_title('Hyper-hydroxymethylated', fontsize=12);

    venn2([set1_down, set2_down], ('Tissue', 'cfDNA'), ax=axs[1]);
    axs[1].set_title('Hypo-hydroxymethylated', fontsize=12);

    return fig

def create_sfigure_3C():
    """Create and save Supplemental Figure 3C"""
    with open(os.path.join(indir, 'SFigure_3C.venn_diagram.pkl'), 'rb') as f:
        venn_diagrams = pickle.load(f)

    fig, axs = plt.subplots(figsize=(8, 6), ncols=3, nrows=2)
    fig.subplots_adjust(wspace=0, hspace=-0.1)
    axs = axs.flatten()
    for i, group in enumerate(['all', 'breast', 'colon', 'lung', 'ovary', 'pancreas']):
        cfdna_size, tissue_size, intersect_size = venn_diagrams[group]
        joint_space = cfdna_size + tissue_size - intersect_size
        M = 189  ## Total number of c6.oncogenic pathways
        table = pd.DataFrame(0, index=[False, True], columns=[False, True])
        table.loc[False, False] = joint_space
        table.loc[False, True] = cfdna_size - intersect_size
        table.loc[True, False] = tissue_size - intersect_size
        table.loc[True, True] = intersect_size
        rv, pval = prob_intersect_hypergeom(M, table)
        
        venn2((venn_diagrams[group]), set_labels=('tissue', 'cfDNA'), ax=axs[i]);
        axs[i].set_title('{:}\np: {:.4f}'.format(group, pval), y=0.9);
    return fig
    
def create_sfigure_5A():
    """Create and save Supplemental Figure 5A"""
    fn = './data/SFigure_5A.roc_curve.TOTO_model_tumor_tissue.tsv.gz'
    roc_curves_to_plot = pd.read_csv(fn, sep='\t', index_col=0)
    roc_curves_to_plot_cancer = roc_curves_to_plot[roc_curves_to_plot['cohort'] == 'Cancer']
    roc_curves_to_plot_normal = roc_curves_to_plot[roc_curves_to_plot['cohort'] == 'Normal']
    palette = roc_curves_to_plot.set_index('label')['color'].to_dict()

    fig, axs = plt.subplots(figsize=(9, 4), ncols=2, sharey=True, sharex=True)
    fig.subplots_adjust(wspace=0.02);
    sns.lineplot(data=roc_curves_to_plot_cancer, x='1 - specificity', y='sensitivity', hue='label', palette=palette, ax=axs[0]);
    axs[0].legend(frameon=True, fontsize=8, ncols=2);
    sns.lineplot(data=roc_curves_to_plot_normal, x='1 - specificity', y='sensitivity', hue='label', palette=palette, ax=axs[1]);
    axs[1].legend(frameon=True, fontsize=8, ncols=2);
    axs[0].set_title('Tumor Tissues', fontsize=10);
    axs[1].set_title('Normal Tissues', fontsize=10);
    add_identity(axs[0], xlim=(0, 1), ylim=(0, 1));
    add_identity(axs[1], xlim=(0, 1), ylim=(0, 1));
    return fig

def create_sfigure_5B():
    """Create and save Supplemental Figure 5B"""
    fn = './data/SFigure_5B.confusion_matrix_top_two.tissue_toto_on_tumor_tissues.tsv.gz'
    df_top_two = pd.read_csv(fn, sep='\t', index_col=0)

    sensitivity = np.diag(df_top_two.values) / df_top_two.sum(0)
    xticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(sensitivity.index, sensitivity.values)]
    specificity = calculate_specificity(df_top_two)
    yticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(specificity.index, specificity.values)][::-1]
    acc = np.diag(df_top_two).sum() / np.sum(df_top_two.values) * 100

    fig, ax = plt.subplots(figsize=(5, 4.5))
    gfg = sns.heatmap(df_top_two.iloc[::-1], annot=True, cmap='Blues', fmt='d', cbar=False, ax=ax, annot_kws={'fontsize': 12}, vmin=0, vmax=50);
    ax.set_title('Tissue TOTO', fontsize=12);
    ax.set_xlabel('True Class (sensitivity)', fontsize=12, labelpad=8);
    ax.set_ylabel('Predicted Class (specificity)', fontsize=12, labelpad=8);
    ax.set_xticklabels(xticks, ha='center', fontsize=10);
    ax.set_yticklabels(yticks, ha='center', fontsize=10, x=-0.12);
    return fig

def create_sfigure_5D():
    """Create and save Supplemental Figure 5D"""
    fn = './data/SFigure_5D.confusion_matrix_top_one.cfdna_toto.tsv.gz'
    df_top_one = pd.read_csv(fn, sep='\t', index_col=0)

    sensitivity = np.diag(df_top_one.values) / df_top_one.sum(0)
    xticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(sensitivity.index, sensitivity.values)]
    specificity = calculate_specificity(df_top_one)
    yticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(specificity.index, specificity.values)][::-1]

    fig, ax = plt.subplots(figsize=(5, 4.5))
    gfg = sns.heatmap(df_top_one.iloc[::-1], annot=True, cmap='Blues', fmt='d', cbar=False, ax=ax, annot_kws={'fontsize': 12}, vmin=0, vmax=50);
    ax.set_xlabel('True Class (sensitivity)', fontsize=12, labelpad=8);
    ax.set_ylabel('Predicted Class (specificity)', fontsize=12, labelpad=8);
    ax.set_xticklabels(xticks, ha='center', fontsize=10);
    ax.set_yticklabels(yticks, ha='center', fontsize=10, x=-0.03);
    return fig

def create_sfigure_5E():
    df = pd.read_csv(os.path.join(indir, 'SFigure_5E_data.csv'), index_col=0)
    TPs = df[df['TPs']]
    df_accuracy = df.groupby('stage_current')['TPs'].sum() / df.groupby('stage_current').size()

    # Create contingency table for chi-square test
    contingency_table = pd.DataFrame({
        'TPs': df.groupby('stage_current')['TPs'].sum(),
        'non_TPs': df.groupby('stage_current').size() - df.groupby('stage_current')['TPs'].sum()
    })
   
    # Perform chi-square test
    chi2, p_value, dof, expected = sp.stats.chi2_contingency(contingency_table, correction=False)
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(6, 3))
    # Create the barplot
    ax.bar(df_accuracy.index, df_accuracy.values, color='steelblue', zorder=10)
    # Add grid lines
    ax.grid(axis='y', linestyle='-', alpha=0.2)
    # Add labels and title
    ax.set_xlabel('Stage', fontsize=12)
    ax.set_ylabel('Accuracy', fontsize=12)
    ax.text(0.5, 0.9, f'chi-square test p={p_value:.3f}',
            fontsize=12,
            transform=ax.transAxes,
            ha='center',
            va='center')
    # Set y-axis limits
    ax.set_ylim(-0.025, 1.025);
    
    return fig

def create_sfigure_5F():
    """Create and save Supplemental Figure 5F"""
    fn = './data/SFigure_5F.confusion_matrix_top_one.cfdna_toto_on_tumor_tissues.tsv.gz'
    df_top_one = pd.read_csv(fn, sep='\t', index_col=0)
    
    sensitivity = np.diag(df_top_one.values) / df_top_one.sum(0)
    xticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(sensitivity.index, sensitivity.values)]
    specificity = calculate_specificity(df_top_one)
    yticks = ['{:}\n({:.1f}%)'.format(x, y*100) for x, y in zip(specificity.index, specificity.values)][::-1]
    acc = np.diag(df_top_one).sum() / np.sum(df_top_one.values) * 100

    fig, ax = plt.subplots(figsize=(5, 4.5))
    gfg = sns.heatmap(df_top_one.iloc[::-1], annot=True, cmap='Blues', fmt='d', cbar=False, ax=ax, annot_kws={'fontsize': 12}, vmin=0, vmax=50);
    ax.set_xlabel('True Class (sensitivity)', fontsize=12, labelpad=8);
    ax.set_ylabel('Predicted Class (specificity)', fontsize=12, labelpad=8);
    ax.set_xticklabels(xticks, ha='center', fontsize=10);
    ax.set_yticklabels(yticks, ha='center', fontsize=10, x=-0.03);
    return fig

#%%```python

if __name__ == "__main__":
    # Create Supplemental Figure 1A
    sfig_1A = create_sfigure_1A() 
    sfig_1A.savefig(fd+'SFigure1A.pdf', bbox_inches='tight')

    # Create Supplemental Figure 2A
    sfig_2A = create_sfigure_2A()
    sfig_2A.savefig(fd+'SFigure2A.pdf', bbox_inches='tight')
    
    sfig_3A = create_sfigure_3A()
    sfig_3A.savefig(fd+'SFigure3A.pdf', bbox_inches='tight')
    
    sfig_3B = create_sfigure_3B()
    sfig_3B.savefig(fd+'SFigure3B.pdf', bbox_inches='tight')
    
    sfig_3C = create_sfigure_3C()
    sfig_3C.savefig(fd+'SFigure3C.pdf', bbox_inches='tight')

    sfig_5A = create_sfigure_5A()
    sfig_5A.savefig(fd+'SFigure5A.pdf', bbox_inches='tight')
    
    sfig_5B = create_sfigure_5B()
    sfig_5B.savefig(fd+'SFigure5B.pdf', bbox_inches='tight')
    
    sfig_5D = create_sfigure_5D()
    sfig_5D.savefig(fd+'SFigure5D.pdf', bbox_inches='tight')
    
    sfig_5E = create_sfigure_5E()
    sfig_5E.savefig(fd+'SFigure5E.pdf', bbox_inches='tight')

    sfig_5F = create_sfigure_5F()
    sfig_5F.savefig(fd+'SFigure5F.pdf', bbox_inches='tight')