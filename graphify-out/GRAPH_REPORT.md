# Graph Report - DeepPBI-KG  (2026-11-06)

## Corpus Check
- Corpus is ~24,371 words - fits in a single context window. You may not need a graph.

## Summary
- 250 nodes · 474 edges · 13 communities (10 shown, 3 thin omitted)
- Extraction: 97% EXTRACTED · 3% INFERRED · 0% AMBIGUOUS · INFERRED: 16 edges (avg confidence: 0.86)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Main prediction pipeline
- Key-gene selection helpers
- KGDeepPBI duplicate pipeline
- Figures and model framework
- Model training and prediction
- Key-gene feature generation
- WGS feature generation
- Setup and BLAST workflow
- WGS vs KG similarity analysis
- GenBank preprocessing
- UniProt to GO mapping
- CDD BLAST script
- Prokka batch script

## God Nodes (most connected - your core abstractions)
1. `protein_process()` - 9 edges
2. `protein_process()` - 8 edges
3. `protein_process()` - 7 edges
4. `protein_process()` - 7 edges
5. `Feature Similarity Report (WGS vs Key Genes)` - 7 edges
6. `compare_features()` - 7 edges
7. `dna_process()` - 6 edges
8. `Dataset` - 6 edges
9. `compare_features()` - 6 edges
10. `model directory listing` - 6 edges

## Surprising Connections (you probably didn't know these)
- `visualize_features.py (planned script)` --conceptually_related_to--> `compare_features()`  [INFERRED]
  .agent/rules/AGENT.md → feature_similarity_report.md
- `WGS vs Key Genes Visualization Implementation Plan` --rationale_for--> `Feature Similarity Report (WGS vs Key Genes)`  [INFERRED]
  .agent/rules/AGENT.md → feature_similarity_report.md
- `Example dataset (9 phages, 8 bacteria)` --shares_data_with--> `template/feature_file intermediate features`  [EXTRACTED]
  .agent/rules/AGENT.md → ReadMe.md
- `Bacterium WGS vs Key Genes Feature Comparison` --semantically_similar_to--> `Phage WGS vs Key Genes Feature Comparison`  [INFERRED] [semantically similar]
  figure/bacterium_wgs_vs_kg_comparison.png → figure/phage_wgs_vs_kg_comparison.png
- `phage_key_gene_0.0001.csv and host_key_gene_0.00004.csv` --semantically_similar_to--> `Top m phage and top n host key gene lists`  [INFERRED] [semantically similar]
  figure/model_parameter.png → figure/model framework.png

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **DeepPBI-KG prediction pipeline** — prokka_blast, deeppbi_kg, readme_result_output, readme_model_parameters [EXTRACTED 1.00]
- **WGS vs KG comparison methods** — feature_similarity_report_cosine_pearson, feature_similarity_report_pca_pc1, feature_similarity_report_joint_pca, feature_similarity_report_heatmap [EXTRACTED 1.00]
- **DeepPBI-KG run workflow (command, template, model, output)** — figure_run_example_command, figure_intermediate_result_example_template_dir, figure_model_parameter_model_dir, figure_output_result_example_output_dir [INFERRED 0.85]
- **Key gene to DNN prediction pipeline** — figure_model_framework_key_gene_vectors, figure_model_framework_random_forest, figure_model_framework_top_key_gene_lists, figure_model_framework_feature_matrix, figure_model_framework_dnn [EXTRACTED 1.00]

## Communities (13 total, 3 thin omitted)

### Community 0 - "Main prediction pipeline"
Cohesion: 0.07
Nodes (31): cds_dna_protein_generate(), check_blank_data(), Count1(), CTDC(), CTDT(), dna_features(), dna_process(), DNN (+23 more)

### Community 1 - "Key-gene selection helpers"
Cohesion: 0.10
Nodes (13): walkFile(), walkFile(), format_fasta(), walkFile(), get_gene(), walkFile(), distance(), Count1() (+5 more)

### Community 2 - "KGDeepPBI duplicate pipeline"
Cohesion: 0.11
Nodes (18): cds_dna_protein_generate(), check_blank_data(), Count1(), CTDC(), CTDT(), dna_features(), dna_process(), DNN (+10 more)

### Community 3 - "Figures and model framework"
Cohesion: 0.09
Nodes (28): Key genes (KG) representation, Whole genome sequence (WGS) representation, Bacterium cross-organism WGS vs KG cosine similarity heatmap, Bacterium WGS vs Key Genes Feature Comparison, Bacterium PC1 WGS vs KG Pearson r=0.9854, template directory (feature_file, host_dna, host_protein, phage_dna, phage_protein), Deep Neural Network output probability, Phage/host DNA and protein feature encoding (+20 more)

### Community 4 - "Model training and prediction"
Cohesion: 0.12
Nodes (3): DNN, Dataset, DNN

### Community 5 - "Key-gene feature generation"
Cohesion: 0.15
Nodes (11): Count1(), CTDC(), CTDT(), dna_features(), dna_process(), protein_features(), protein_process(), walkFile() (+3 more)

### Community 6 - "WGS feature generation"
Cohesion: 0.16
Nodes (9): Count1(), CTDC(), CTDT(), dna_features(), dna_process(), protein_features(), protein_process(), walkFile() (+1 more)

### Community 7 - "Setup and BLAST workflow"
Cohesion: 0.16
Nodes (13): walkFile(), prokka_blast.sh script, BLAST alignment, Conda environment DeepPBI-KG, DeepPBI-KG, all_host_db BLAST database, all_phage_db BLAST database, Prokka annotation (+5 more)

### Community 8 - "WGS vs KG similarity analysis"
Cohesion: 0.27
Nodes (11): Example dataset (9 phages, 8 bacteria), WGS vs Key Genes Visualization Implementation Plan, visualize_features.py (planned script), compare_features(), Per-organism cosine similarity and Pearson correlation, Cross-organism cosine similarity heatmap, Joint 2D PCA with WGS-KG connections, PCA PC1 linear relationship (+3 more)

### Community 9 - "GenBank preprocessing"
Cohesion: 0.31
Nodes (5): format_fasta(), get_cds_dna(), get_cds_protein(), get_complete(), walkFile()

## Knowledge Gaps
- **11 isolated node(s):** `Conda environment DeepPBI-KG`, `Raw vs scaled (Z-score) features`, `torch 1.9.0+cpu`, `scikit-learn 0.23.2`, `biopython 1.78` (+6 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 73 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **3 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `template/feature_file intermediate features` connect `WGS vs KG similarity analysis` to `Main prediction pipeline`?**
  _High betweenness centrality (0.054) - this node is a cross-community bridge._
- **What connects `Conda environment DeepPBI-KG`, `Raw vs scaled (Z-score) features`, `torch 1.9.0+cpu` to the rest of the system?**
  _11 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Main prediction pipeline` be split into smaller, more focused modules?**
  _Cohesion score 0.06859903381642513 - nodes in this community are weakly interconnected._
- **Why does `compare_features()` connect `WGS vs KG similarity analysis` to `Main prediction pipeline`?**
  _High betweenness centrality (0.037) - this node is a cross-community bridge._
- **Should `Key-gene selection helpers` be split into smaller, more focused modules?**
  _Cohesion score 0.09716599190283401 - nodes in this community are weakly interconnected._
- **Why does `Prediction result (result.csv / result.html)` connect `Setup and BLAST workflow` to `Main prediction pipeline`?**
  _High betweenness centrality (0.035) - this node is a cross-community bridge._
- **Should `KGDeepPBI duplicate pipeline` be split into smaller, more focused modules?**
  _Cohesion score 0.1103448275862069 - nodes in this community are weakly interconnected._