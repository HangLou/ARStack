from scipy import stats
import pandas as pd
import seaborn as sns
from statsmodels.stats.weightstats import ztest
import matplotlib.pyplot as plt
dataset_rmsd1 = pd.read_csv('Dataset1_stack_alpha_baker_boxplot.csv')


def significance_test(df, data_name='Dataset 1'):
    ks_stack_alpha = stats.ks_2samp(df['stack'].values, df['alpha'].values)
    t_stack_alpha = stats.ttest_ind(df['stack'].values, df['alpha'].values)
    ks_stack_baker = stats.ks_2samp(df['stack'].values, df['baker'].values)
    t_stack_baker = stats.ttest_ind(df['stack'].values, df['baker'].values)
    pivot = df.iloc[:, 1:4].stack().reset_index()
    pivot = pivot.rename(columns={'level_1': 'model', 0: 'rmsd'})
    sns.set()

    plt.figure()
    rmsd_dist = sns.displot(
        pivot, x='rmsd', hue='model', kind='kde', alpha=0.1)

    #rmsd_dist.set_title('RMSD distribution comparison on '+data_name)
    #fig = rmsd_dist.get_figure()
    plt.title('RMSD distribution comparison on '+data_name)
    plt.tight_layout()
    plt.savefig('rmsd dist_plot'+data_name+".png")
    print('p-value for  betweent stack model and alphafold=',
          ks_stack_alpha[1])
    print('p-value for T test betweent stack model and alphafold=',
          t_stack_alpha[1])
    print('p-value for Kolmogorov-Smirnov test betweent stack model and baker=',
          ks_stack_baker[1])
    print('p-value for T test between stack model and baker=',
          t_stack_baker[1])
    df['diff'] = df.iloc[:, 2]-df.iloc[:, 1]
    print(df)
    df = df.rename(columns={'diff': '\u0394RMSD'})
    plt.figure()
    ax = sns.kdeplot(df.iloc[:, -1])
    line = ax.get_lines()[-1]
    x, y = line.get_data()
    mask = x >= 0
    x, y = x[mask], y[mask]
    # print(mask)
    ax.fill_between(x, y1=y, alpha=0.3, facecolor='green')

    line = ax.get_lines()[-1]
    x, y = line.get_data()
    mask = x <= 0.0005
    x, y = x[mask], y[mask]
    # print(x)
    ax.fill_between(x, y1=y, alpha=0.3, facecolor='red')
    plt.title('\u0394RMSD distribution plot on '+data_name)
    #drmsd_dist.set_title('\u0394RMSD distribution plot on '+data_name)
    #fig = drmsd_dist.get_figure()
    plt.tight_layout()
    plt.savefig('drmsd dist_plot'+data_name+".png")

    plt.figure()
    drmsd_boxplot = sns.boxplot(df.iloc[:, -1])
    #drmsd_boxplot.set_title('\u0394RMSD boxplot on '+data_name)
    #fig = drmsd_boxplot.get_figure()

    plt.title('\u0394RMSD boxplot on '+data_name)
    plt.tight_layout()
    plt.savefig('drmsd boxplot_plot'+data_name+".png")

    z_test = ztest(df.iloc[:, -1].values, alternative='smaller')
    print('p-value for z test between on the dRMSD between stack model and alpha=',
          z_test[1])


def plot():
    sns.set()
    dataset_rmsd1 = pd.read_csv('Dataset1_stack_alpha_baker_boxplot.csv')
    dataset_rmsd2 = pd.read_csv('Dataset2_stack_alpha_baker_boxplot.csv')
    dataset_rmsd3 = pd.read_csv('Dataset3_stack_alpha_baker_boxplot.csv')
    dataset_rmsd1['dataset'] = 1
    dataset_rmsd2['dataset'] = 2
    dataset_rmsd3['dataset'] = 3
    MW1 = pd.read_csv('Dataset_1_MW.csv', names=['MW'])
    MW2 = pd.read_csv('Dataset_2_MW.csv', names=['MW'])
    MW3 = pd.read_csv('Dataset_3_MW.csv', names=['MW'])
    MW = pd.concat([MW1, MW2, MW3])
    MW = MW.reset_index()
    MW = MW.rename(columns={'index': 'protein'})
    df = pd.concat([dataset_rmsd1, dataset_rmsd2, dataset_rmsd3])
    df['diff'] = df.iloc[:, 2]-df.iloc[:, 1]

    df = df.rename(columns={'diff': '\u0394RMSD'})
    df = df.rename(columns={'Unnamed: 0': 'protein'})
    df = pd.merge(df, MW, how='left',
                  on='protein').drop_duplicates(keep='first')

    plt.figure(figsize=(8, 6))
    ax = sns.violinplot(data=df, x='dataset',
                        y='\u0394RMSD', linewidth=2.5)
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("drmsd_violine_plot.png")

    df['temp'] = df.MW <= 10000
    df['Molecular Weight (W)'] = r'$W \leq 10000$'
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$W \leq 10000$'

    df['temp'] = (df.MW > 10000) & (df.MW <= 20000)
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$10000 < W \leq 20000$'

    df['temp'] = df.MW > 20000
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$20000 < W$'

    plt.figure(figsize=(8, 6))

    ax = sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
                        hue='Molecular Weight (W)', linewidth=2.5)
    #handles, labels = ax.get_legend_handles_labels()

    # ax.legend(handles, [r'$W<=10000$', r'10000<W<=20000',
    #         r'$20000<W$'], loc='upper right')
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("drmsd_violine_plot_MW.png")

    #df['Alphafold RMSD > 5'] = df.alpha > 5
    # sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
    #              hue='Alphafold RMSD > 5')
    plt.figure(figsize=(8, 6))
    df['Alphafold'] = df.alpha > 1
    df.Alphafold.loc[df.Alphafold == True] = r'$RMSD > 1$'
    df.Alphafold.loc[df.Alphafold == False] = r'$RMSD \leq 1$'
    ax = sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
                        hue='Alphafold', linewidth=2.5)
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("drmsd_violine_plot_cohorts.png")


if __name__ == '__main__':
    dataset_rmsd1 = pd.read_csv('Dataset1_stack_alpha_baker_boxplot.csv')
    dataset_rmsd2 = pd.read_csv('Dataset2_stack_alpha_baker_boxplot.csv')
    dataset_rmsd3 = pd.read_csv('Dataset3_stack_alpha_baker_boxplot.csv')
    significance_test(dataset_rmsd1)
