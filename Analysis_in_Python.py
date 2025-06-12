import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import mannwhitneyu

# Constants
GTF_PATH = "d:/Learning/Drosofila/genomic.gtf"
COUNTS_PATH = "d:/Learning/Drosofila/filtered_counts.txt"
GENE_INFO_PATH = "d:/Learning/Drosofila/List_of_genes_with_gene_id.csv"
GROUPS = {
    'male': ['mOregon_1', 'mOregon_2', 'mCHD1_1', 'mCHD1_2'],
    'female': ['fOregon_1', 'fOregon_2', 'fCHD1_1', 'fCHD1_2']
}
ALL_SAMPLES = GROUPS['female'] + GROUPS['male']

# Functions
def parse_gtf_lengths(gtf_file):
    gene_lengths = defaultdict(int)
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "exon":
                continue
            start, end = int(fields[3]), int(fields[4])
            info = fields[8]
            gene_id = [s.split('"')[1] for s in info.split(';') if 'gene_id' in s][0]
            gene_lengths[gene_id] += (end - start + 1)
    return pd.DataFrame(gene_lengths.items(), columns=['gene_id', 'length'])

def counts_to_tpm(df, sample_cols):
    tpm_df = df[['gene_id']].copy()
    for sample in sample_cols:
        rpk = df[sample] / df['length']
        tpm_df[sample] = rpk / (rpk.sum() / 1e6)
    return tpm_df

def process_sex(tpm, sex, hk_genes):
    samples = GROUPS[sex]
    cond_samples = [s for s in samples if 'CHD1' in s]
    control_samples = [s for s in samples if 'Oregon' in s]

    # Relaxed filter: TPM > 1
    tpm_sex = tpm[(tpm[samples] > 1).all(axis=1)].copy()
    tpm_sex[f'{sex}_CHD1_mean'] = tpm_sex[cond_samples].mean(axis=1)
    tpm_sex[f'{sex}_Oregon_mean'] = tpm_sex[control_samples].mean(axis=1)

    ratio = tpm_sex[f'{sex}_CHD1_mean'] / tpm_sex[f'{sex}_Oregon_mean']
    ratio.replace(0, np.nan, inplace=True)
    tpm_sex['log2FC'] = np.log2(ratio)
    tpm_sex.dropna(subset=['log2FC'], inplace=True)

    # Centering log2FC
    tpm_sex['log2FC_centered'] = tpm_sex['log2FC'] - tpm_sex['log2FC'].mean()
    tpm_sex = tpm_sex[(tpm_sex['log2FC_centered'] >= -1) & (tpm_sex['log2FC_centered'] <= 1)]

    annotated = pd.merge(tpm_sex, hk_genes[['gene_id', 'Gene name', 'Chr']], on='gene_id', how='left')
    annotated['Sex'] = sex.capitalize()
    return annotated

def plot_box_density(data, sex, suffix):
    chrX = data[data['Chr'] == 'X']['log2FC_centered']
    autosomes = data[data['Chr'] != 'X']['log2FC_centered']

    # Boxplot
    plt.figure(figsize=(8, 6))
    plt.boxplot([chrX.dropna(), autosomes.dropna()], tick_labels=['chrX', 'Autosomes'])
    plt.axhline(0, color='gray', linestyle='--')
    plt.title(f"Dosage Compensation (Housekeeping Genes) - {sex}")
    plt.ylabel("Centered log2FC")
    plt.tight_layout()
    plt.savefig(f"boxplot_{suffix}.png")
    plt.close()

    # KDE plot
    plt.figure(figsize=(8, 6))
    sns.kdeplot(chrX, label='chrX', color='orange', linewidth=2)
    sns.kdeplot(autosomes, label='Autosomes', color='green', linewidth=2)
    plt.axvline(chrX.median(), color='orange', linestyle='--', label=f'Median chrX: {chrX.median():.2f}')
    plt.axvline(autosomes.median(), color='green', linestyle='--', label=f'Median Autosomes: {autosomes.median():.2f}')
    plt.axvline(0, color='gray', linestyle=':')
    plt.legend()
    plt.title(f"Distribution of Centered log2FC - {sex}")
    plt.xlabel("Centered log2FC")
    plt.tight_layout()
    plt.savefig(f"density_{suffix}.png")
    plt.close()

# Load count data
with open(COUNTS_PATH) as f:
    data_lines = [line.strip() for line in f if not line.startswith("#") and not line.startswith("Geneid")]
counts = pd.DataFrame([line.split() for line in data_lines], columns=['gene_id'] + ALL_SAMPLES)
counts[ALL_SAMPLES] = counts[ALL_SAMPLES].astype(int)

# Load housekeeping gene list
hk_genes = pd.read_csv(GENE_INFO_PATH, sep=";")
counts = counts[counts['gene_id'].isin(hk_genes['gene_id'])]

# Parse gene lengths
lengths = parse_gtf_lengths(GTF_PATH)
df = pd.merge(counts, lengths, on="gene_id")
tpm = counts_to_tpm(df, ALL_SAMPLES)

# Analysis and visualization
combined_results = []
for sex in ['male', 'female']:
    annotated = process_sex(tpm, sex, hk_genes)

    # Mann–Whitney U test
    chrX_vals = annotated[annotated['Chr'] == 'X']['log2FC_centered'].dropna()
    autosome_vals = annotated[annotated['Chr'] != 'X']['log2FC_centered'].dropna()
    stat, pval = mannwhitneyu(chrX_vals, autosome_vals, alternative='two-sided')

    print(f"{sex.capitalize()} - Mann–Whitney U test p-value: {pval:.4e}")

    plot_box_density(annotated, f"{sex.capitalize()} (p = {pval:.2e})", f"{sex}_HK_chrX_vs_auto")

    filename = f"filtered_housekeeping_genes_log2FC_centered_{sex}_sorted.csv"
    annotated[['Gene name', 'log2FC_centered']].rename(
        columns={'log2FC_centered': f'log2FC_centered_{sex}'}
    ).sort_values(by=f'log2FC_centered_{sex}', ascending=False).to_csv(filename, index=False)

    combined_results.append(annotated)

# Combine male and female results
combined = pd.concat(combined_results)
combined['ChrType'] = combined['Chr'].apply(lambda x: 'chrX' if x == 'X' else 'Autosomes')

# KDE plots by chromosome type and sex
for chr_type in ['X', '!= X']:
    chr_label = 'chrX' if chr_type == 'X' else 'Autosomes'
    chr_data = combined.query(f'Chr {"==" if chr_type == "X" else "!="} "X"')
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=chr_data, x='log2FC_centered', hue='Sex', common_norm=False, linewidth=2)
    for sex in ['Female', 'Male']:
        med = chr_data[chr_data['Sex'] == sex]['log2FC_centered'].median()
        color = 'orange' if sex == 'Female' else 'blue'
        plt.axvline(med, color=color, linestyle='--', label=f'Median {sex} {chr_label}: {med:.2f}')
    plt.axvline(0, color='gray', linestyle=':')
    plt.title(f"Distribution of Centered log2FC on {chr_label} by Sex")
    plt.xlabel("Centered log2FC")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"kde_log2FC_centered_{chr_label}_by_sex.png")
    plt.close()

# Final boxplot by sex and chromosome type
plt.figure(figsize=(10, 6))
sns.boxplot(data=combined, x='ChrType', y='log2FC_centered', hue='Sex')
plt.axhline(0, color='gray', linestyle='--')
plt.title("Boxplot of Centered log2FC by Sex and Chromosome Type")
plt.xlabel("Chromosome Type")
plt.ylabel("Centered log2FC")
plt.legend(title="Sex")
plt.tight_layout()
plt.savefig("boxplot_log2FC_centered_by_sex_and_chrType.png")
plt.close()
