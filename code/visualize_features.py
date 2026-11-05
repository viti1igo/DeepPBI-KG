import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr

# Set plot style for premium aesthetics
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.titlesize': 18
})

# Palettes
C_WGS = "#1f77b4" # Blue
C_KG = "#ff7f0e"  # Orange

def compare_features(wgs_path, kg_path, scaler_path, species_name, output_prefix):
    print(f"--- Analyzing {species_name} ---")
    
    # Load features
    df_wgs = pd.read_csv(wgs_path, index_col=0)
    df_kg = pd.read_csv(kg_path, index_col=0)
    
    # Ensure indices align
    common_idx = df_wgs.index.intersection(df_kg.index)
    df_wgs = df_wgs.loc[common_idx]
    df_kg = df_kg.loc[common_idx]
    
    names = common_idx.tolist()
    print(f"Loaded {len(names)} samples. Feature dimensions: {df_wgs.shape[1]}")
    
    # Load scaler and extract species-specific parameters
    # phage (1338) + host (1338) = 2676 features
    with open(scaler_path, 'rb') as f:
        scaler = pickle.load(f)
    
    mean_val = scaler.mean_
    scale_val = scaler.scale_
    
    if species_name.lower() == 'phage':
        mean_sp = mean_val[:1338]
        scale_sp = scale_val[:1338]
    else:
        mean_sp = mean_val[1338:]
        scale_sp = scale_val[1338:]
        
    # Scale features
    feats_wgs_scaled = (df_wgs.values - mean_sp) / scale_sp
    feats_kg_scaled = (df_kg.values - mean_sp) / scale_sp
    
    # Compute Cosine Similarity and Pearson Correlation for each sample
    cosine_sims_raw = []
    cosine_sims_scaled = []
    pearson_corrs_raw = []
    pearson_corrs_scaled = []
    
    for i in range(len(names)):
        v_wgs_raw = df_wgs.values[i]
        v_kg_raw = df_kg.values[i]
        v_wgs_sc = feats_wgs_scaled[i]
        v_kg_sc = feats_kg_scaled[i]
        
        # Raw Cosine Similarity
        cos_raw = np.dot(v_wgs_raw, v_kg_raw) / (np.linalg.norm(v_wgs_raw) * np.linalg.norm(v_kg_raw) + 1e-9)
        cosine_sims_raw.append(cos_raw)
        
        # Scaled Cosine Similarity
        cos_sc = np.dot(v_wgs_sc, v_kg_sc) / (np.linalg.norm(v_wgs_sc) * np.linalg.norm(v_kg_sc) + 1e-9)
        cosine_sims_scaled.append(cos_sc)
        
        # Pearson Correlation
        p_raw, _ = pearsonr(v_wgs_raw, v_kg_raw)
        pearson_corrs_raw.append(p_raw)
        
        p_sc, _ = pearsonr(v_wgs_sc, v_kg_sc)
        pearson_corrs_scaled.append(p_sc)
        
    metrics_df = pd.DataFrame({
        'Organism': names,
        'Cosine_Similarity_Raw': cosine_sims_raw,
        'Cosine_Similarity_Scaled': cosine_sims_scaled,
        'Pearson_Correlation_Raw': pearson_corrs_raw,
        'Pearson_Correlation_Scaled': pearson_corrs_scaled
    })
    print("\nSimilarity Metrics summary:")
    print(metrics_df.describe().to_string())
    
    # Dimensionality Reduction (PCA on joint WGS & KG representations)
    combined_scaled = np.vstack([feats_wgs_scaled, feats_kg_scaled])
    pca = PCA(n_components=2)
    pca.fit(combined_scaled)
    
    coords_wgs = pca.transform(feats_wgs_scaled)
    coords_kg = pca.transform(feats_kg_scaled)
    
    # Create visualization plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.suptitle(f"{species_name} WGS vs. Key Genes (KG) Feature Comparison", y=0.98)
    
    # Subplot A: PC1 WGS vs PC1 KG (Linear relationship test)
    ax_lin = axes[0, 0]
    pc1_wgs = coords_wgs[:, 0]
    pc1_kg = coords_kg[:, 0]
    r_val, p_val = pearsonr(pc1_wgs, pc1_kg)
    
    sns.regplot(x=pc1_wgs, y=pc1_kg, ax=ax_lin, color="purple",
                scatter_kws={"s": 80, "alpha": 0.8, "edgecolor": "k"})
    ax_lin.set_xlabel("WGS Representation (PC1)")
    ax_lin.set_ylabel("Key Genes Representation (PC1)")
    ax_lin.set_title(f"Linear Relationship of PC1\n(Pearson r = {r_val:.4f}, p = {p_val:.1e})")
    
    # Subplot B: 2D PCA Space with connections
    ax_pca = axes[0, 1]
    # Plot WGS points
    ax_pca.scatter(coords_wgs[:, 0], coords_wgs[:, 1], color=C_WGS, label="WGS (Whole Genome)", s=90, edgecolors='k', zorder=3)
    # Plot KG points
    ax_pca.scatter(coords_kg[:, 0], coords_kg[:, 1], color=C_KG, label="KG (Key Genes)", s=90, edgecolors='k', marker='D', zorder=3)
    
    # Draw connections
    for i in range(len(names)):
        ax_pca.plot([coords_wgs[i, 0], coords_kg[i, 0]], [coords_wgs[i, 1], coords_kg[i, 1]], 
                    color="gray", linestyle="--", alpha=0.6, zorder=2)
        # Add labels
        ax_pca.text(coords_wgs[i, 0] + 0.1, coords_wgs[i, 1] + 0.1, names[i], fontsize=9, alpha=0.8)
        
    ax_pca.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    ax_pca.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    ax_pca.set_title("Joint 2D PCA Space with WGS-KG Connections")
    ax_pca.legend(frameon=True)
    
    # Subplot C: Cosine Similarity and Correlation bar plot (scaled features)
    ax_bar = axes[1, 0]
    x_indices = np.arange(len(names))
    width = 0.35
    
    ax_bar.bar(x_indices - width/2, cosine_sims_scaled, width, label='Cosine Sim (Scaled)', color='#2ca02c', edgecolor='black', alpha=0.8)
    ax_bar.bar(x_indices + width/2, pearson_corrs_scaled, width, label='Pearson Corr (Scaled)', color='#d62728', edgecolor='black', alpha=0.8)
    
    ax_bar.set_xticks(x_indices)
    ax_bar.set_xticklabels(names, rotation=45, ha='right')
    ax_bar.set_title("WGS vs. KG Similarity per Organism")
    ax_bar.set_ylabel("Similarity Metric Score")
    ax_bar.set_ylim(0.0, 1.05) # adjust y-limits to show all values from 0 to 1
    ax_bar.legend(loc='lower left', frameon=True)
    
    # Subplot D: Cosine Similarity Heatmap between all WGS and KG samples
    ax_hm = axes[1, 1]
    # Compute cross cosine similarity matrix
    similarity_matrix = np.zeros((len(names), len(names)))
    for idx_w, w_vec in enumerate(feats_wgs_scaled):
        for idx_k, k_vec in enumerate(feats_kg_scaled):
            sim = np.dot(w_vec, k_vec) / (np.linalg.norm(w_vec) * np.linalg.norm(k_vec) + 1e-9)
            similarity_matrix[idx_w, idx_k] = sim
            
    sim_df = pd.DataFrame(similarity_matrix, index=names, columns=names)
    sns.heatmap(sim_df, annot=True, fmt=".3f", cmap="YlGnBu", ax=ax_hm, cbar_kws={'label': 'Cosine Similarity'})
    ax_hm.set_title("Cross-Organism WGS vs. KG Cosine Similarity Heatmap")
    ax_hm.set_xlabel("Key Genes (KG) Representation")
    ax_hm.set_ylabel("Whole Genome (WGS) Representation")
    plt.setp(ax_hm.get_xticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    plot_path = f"figure/{output_prefix}_comparison.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {plot_path}")
    
    return metrics_df, r_val

if __name__ == '__main__':
    os.makedirs('figure', exist_ok=True)
    
    # Phage 
    phage_metrics, phage_r = compare_features(
        wgs_path='template/feature_file/phage_wgs_feature.csv',
        kg_path='template/feature_file/phage_key_gene_feature.csv',
        scaler_path='model/wgs_scaler.pkl', # DeepPBI-KG uses wgs_scaler.pkl or key_gene_scaler.pkl
        species_name='Phage',
        output_prefix='phage_wgs_vs_kg'
    )
    
    # Bacterium 
    bacterium_metrics, bacterium_r = compare_features(
        wgs_path='template/feature_file/bacterium_wgs_feature.csv',
        kg_path='template/feature_file/bacterium_key_gene_feature.csv',
        scaler_path='model/wgs_scaler.pkl', # scaler is trained on combined phage+host vector
        species_name='Bacterium',
        output_prefix='bacterium_wgs_vs_kg'
    )
    