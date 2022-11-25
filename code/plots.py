from scipy import stats
import pandas as pd
import seaborn as sns
dataset_rmsd1 = pd.read_csv('Dataset1_stack_alpha_baker_boxplot.csv')


def significance_test(df):
    ks_stack_alpha = stats.ks_2samp(df['stack'].values, df['alpha'].values)
    t_stack_alpha = stats.ttest_ind(df['stack'].values, df['alpha'].values)
    ks_stack_baker = stats.ks_2samp(df['stack'].values, df['baker'].values)
    t_stack_baker = stats.ttest_ind(df['stack'].values, df['baker'].values)
    pivot = df.iloc[:, 1:4].stack().reset_index()
    pivot = pivot.rename(columns={'level_1': 'model', 0: 'rmsd'})
    sns.displot(pivot, x='rmsd', hue='model', kind='kde', alpha=0.1)

    print('p-value for Kolmogorov-Smirnov test betweent stack model and alphafold=',
          ks_stack_alpha[1])
    print('p-value for T test betweent stack model and alphafold=',
          t_stack_alpha[1])
    print('p-value for Kolmogorov-Smirnov test betweent stack model and baker=',
          ks_stack_baker[1])
    print('p-value for T test betweent stack model and baker=',
          t_stack_baker[1])


if __name__ == '__main__':
    dataset_rmsd1 = pd.read_csv('Dataset1_stack_alpha_baker_boxplot.csv')
    dataset_rmsd2 = pd.read_csv('Dataset2_stack_alpha_baker_boxplot.csv')
    dataset_rmsd3 = pd.read_csv('Dataset3_stack_alpha_baker_boxplot.csv')
    significance_test(dataset_rmsd1)
