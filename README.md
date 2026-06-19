# CEA_PLMs
We developed the Coupled Epistatic Algorithm (CEA) and integrated it with PLMs to uncover functional and adaptive epistatic mutations in IAV HA. It fully covers known binding sites and reveals their interactions with scaffold regions, establishing an integrated framework for predicting host adaptation and pandemic risk.
This package provides an pipeline of data parsing, coupled epistatic algorithm to select coupled epistatic site pairs and pairwise mutations, ESM-2 fine-tuning, vBERT-based classifier establishment, re-selection based on mutation effect and adaptation shift and benchmarking of the integrated workflow.

Due to GitHub's limitations on the number and size of uploaded files, please download the complete dataset and code from Zenodo (https://doi.org/10.5281/zenodo.18437264). Only the main code is uploaded to GitHub. To use the code, please download the data folder and the subfolders and .csv files in other folders from Zenodo first. These files should be placed in the corresponding folder before running the GitHub code. 

## Hardware/Computational resources
Typical install time: The software is provided as a ready-to-run package and requires no compilation. Installation primarily involves setting up the Python environment and installing the required dependencies. On a "normal" desktop computer, this process typically takes 5 to 10 minutes, depending on the internet connection speed.

Expected run time for demo: Executing the analysis on the HA dataset using the pre-trained model takes approximately 30 minutes on a standard desktop computer. Alternatively, training the model from scratch on the same hardware requires approximately 2 to 3 hours.

Major computational experiments were executed on GPU servers. Server 1: The core hardware configuration features an Intel(R) Xeon(R) Gold 6548Q processor coupled with 512 GB of DIMM Synchronous Registered (Buffered) (5600 MHz). Accelerated training and inference of the deep learning framework were performed utilizing an NVIDIA H100 GPU. The analytical pipeline was deployed within a Linux-based operating system. Server 2: The core hardware configuration features an Intel(R) Xeon(R) Gold 6254 processor coupled with 1280 GB of DIMM DDR4 Synchronous Registered (Buffered) (2666 MHz). Accelerated training and inference of the deep learning framework were performed utilizing an V100 GPU. The analytical pipeline was deployed within a Linux-based operating system. Sequence Alignment was conducted with MAFFT (version 7.511); Sequence Clustering used CD-HIT (version 4.8.1). Protein embeddings and mutational effect predictions were generated and fine-tuned using vBERT (implemented in TensorFlow 2.3.0) and ESM-2 (implemented in PyTorch 2.5.1) architectures. Custom Python (version 3.8) scripts were developed to implement the coupled epistatic algorithm (CEA) to construct the synergistic networks, and to train the host adaptation classifier.  

## Data parsing of the IAV HA sequences and alignment with the reference numbering scheme (reference numbering scheme PMID: 25391151)
1. Amino acid sequence translation and cleaning of IAV HA sequences
```bash
python data_parsing/0translate_nt_to_AA_for_H13579.py
```

2. Aligning the 'data/generate/HA_with_aa_seq_without_wrong_seq_H1.fasta', 'data/generate/HA_with_aa_seq_without_wrong_seq_H3.fasta', and 'data/generate/HA_with_aa_seq_without_wrong_seq_H579.fasta' with mafft and the aligned sequences were saved as 'data/mafft/HA_with_aa_seq_without_wrong_seq_H1_mafft.fasta', 'data/mafft/HA_with_aa_seq_without_wrong_seq_H3_mafft.fasta', and 'data/mafft/HA_with_aa_seq_without_wrong_seq_H579_mafft.fasta'.

3. Data cleaning of the aligned sequences and adding the reference sequence for each serotype
```bash
python data_parsing/1clean_H1.py
python data_parsing/1clean_H3.py
python data_parsing/1clean_H5.py
python data_parsing/1clean_H7.py
python data_parsing/1clean_H9.py
```

4. Alignment to the reference numbering scheme
```bash
python data_parsing/2add_empty_with_H1_PlosOne.py
python data_parsing/2add_empty_with_H3_PlosOne.py
python data_parsing/2add_empty_with_H5_PlosOne.py
python data_parsing/2add_empty_with_H7_PlosOne.py
python data_parsing/2add_empty_with_H9_PlosOne.py
```

5. Creating the fasta file with segmented HA sequences from site 70 to 290 without empty for CD-HIT sampling
```bash
python data_parsing/3create_fasta_without_empty_for_cdhit_sampling_H135.py
```

6. CD-HIT sampling was performed to the segmented H1, H3, and H5 HA sequences with a threshold of 0.995. The input fasta files were 'data/generate/H135_cdhit_sampling/H1_aligned_with_PlosOne_before_sampling.fasta', 'data/generate/H135_cdhit_sampling/H3_aligned_with_PlosOne_before_sampling.fasta', and 'data/generate/H135_cdhit_sampling/H5_aligned_with_PlosOne_before_sampling.fasta'. The sampled sequences were saved as 'data/cdhit/H1_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta', 'data/cdhit/H3_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta'
, and 'data/cdhit/H5_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta'.

7. Merging the CD-HIT-sampled sequence with the H7 and H9 sequences.
```bash
python data_parsing/4extract_cdhit_sampled_seqs.py
```
8. Secondary structure prediction was performed to the H1, H3, H5, H7, and H9 reference HA sequences with AlphaFold2 and the 'ranked_0.pdb' file of each sequence was extracted and saved as 'data/pdb_AlphaFold2/ref_H1_ranked_0.pdb', 'data/pdb_AlphaFold2/ref_H3_ranked_0.pdb', 'data/pdb_AlphaFold2/ref_H5_ranked_0.pdb', 'data/pdb_AlphaFold2/ref_H7.pdb' and 'data/pdb_AlphaFold2/ref_H9.pdb'.

9. Secondary structure of each site was extracted by Pymol, and the txt files were saved as 'data/pdb_AlphaFold2/ref_H1_PlosOne_ref.txt', 'data/pdb_AlphaFold2/ref_H3_PlosOne_ref.txt', 'data/pdb_AlphaFold2/ref_H5_Texas_ref.txt', 'data/pdb_AlphaFold2/ref_H7.txt', and 'data/pdb_AlphaFold2/ref_H9.txt'.

10. Extracting the  secondary structure of each site and junction information
```bash
python data/pdb_AlphaFold2/txt_to_csv_marking_junction_and_marking_RBD.py
```

## Coupled Epistatic Algorithm (CEA) design and benchmarking
1. The 2,000 IAV HA sequences in the test dataset with for vBERT benchmarking (reference: PMID: 41035817) was aligned with the DMS reference HA sequence (reference: PMID: 27271655), and the aligned space were cleaned according to the reference sequence. The parsed sequences for CEA benchmarking were saved as 'benchmarking/HA_test_set_amino_with_H1_HA_DMS_WT_mafft_del_space_based_on_WT.fasta'. The 1,000 SARS-CoV-2 Spike RBD sequences in the test dataset with for vBERT benchmarking (reference: PMID: 41035817) was aligned with the DMS reference Spike RBD sequence (reference: PMID: 36535326), and the aligned space were cleaned according to the reference sequence. The parsed sequences for CEA benchmarking were saved as 'benchmarking/Spike_RBD_for_dms_coevo_matrix.fasta'.


2. CEA calculation based on the HA and Spike RBD sequences for benchmarking
```bash
python benchmarking/1coevo_new_matrix_origin_evaluation_HA_for_DMS_benchmarking.py
python benchmarking/1coevo_new_matrix_origin_evaluation_HA_for_DMS_benchmarking.py
```

3. Site importance calculation based on HA  and Spike RBD DMS experimental data
```bash
python benchmarking/2H1_HA_mutants_preference_score_pos_nag.py
python benchmarking/2S_RBD_DMS_binding_score.py
```

4. Mutation effect calculation based on original ESM-2 650M model of HA and Spike RBD DMS reference sequence
```bash
python benchmarking/DMS_esm2_logist_matrix/HA_DMS_logist_mat_each_site.py
python benchmarking/DMS_esm2_logist_matrix/S_RBD_DMS_logist_mat_each_site.py
```

5. Site importance calculation based on mean field Direct Coupling Analysis (mfDCA) of HA and Spike RBD DMS reference sequence (Please first download the pydca in the 'benchmarking/DCA/' folder)
```bash
python benchmarking/DCA/pydca_2_protein_HA_for_DMS_benchmarking.py
python benchmarking/DCA/pydca_2_protein_S_RBD_for_DMS_benchmarking.py
```
6. Benchmarking of CAE against ESM-2 and DCA on HA and Spike RBD DMS dataset
```bash
python benchmarking/3coevo_dca_benchmarking_with_H1_HA_DMS_new_coevo_with_plot.py
python benchmarking/3coevo_dca_benchmarking_with_H1_HA_DMS_new_coevo_with_plot.py
```

## ESM-2 fine-tuning, benchmarking, and mutation effect calculation of the HA sequences (Please first download the original ESM-2_650M and save in the 'ESM-2_ft' folder)
1. ESM-2 fine-tuning with the HA sequence dataset
```bash
python ESM-2_ft/1esm2_finetune_JSY1_cdhit_freeze32_lr5e-6.py
```

2. ESM-2 fine-tuning with vBERT (reference: PMID: 41035817) pre-training dataset (with multi-viral proteins)
```bash
python ESM-2_ft/2esm2_finetune_JSY1_virus10w_freeze32_lr3e-5.py
```

3. Mutation effect calculation of the DMS reference HA sequence by fine-tuned ESM-2 based on  ESM-2 fine-tuning with the HA sequence dataset
```bash
python ESM-2_ft/3esm2_mutation_effect_masking_JSY_ft1_cdhit999_freeze32_lr5e-6.py
```

4. Mutation effect calculation of the DMS reference HA sequence by fine-tuned ESM-2 based on  ESM-2 fine-tuning with multi-viral proteins
```bash
python ESM-2_ft/4esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
```
5. Fine-tuned ESM-2 benchmarking based on HA DMS experimental dataset
```bash
python benchmarking/4H1_HA_mutants_site_entropy.py
python benchmarking/5DMS_for_esm2_score.py
python benchmarking/6_1DMS_minus_origin_AA_HA_and_spearman_with_esm2_ft_models_HA_compare_with_plot.py
python benchmarking/6_2DMS_minus_origin_AA_HA_and_spearman_with_esm2_ft_models_with_plot.py
```

6. Mutation effect calculation of the H1,H3, H5, H7, and H9 reference HA sequences by fine-tuned ESM-2 based on  ESM-2 fine-tuning with multi-viral proteins (selected based on benchmarking)
```bash
python ESM-2_ft/5CA07_H1_ref_esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
python ESM-2_ft/6PlosOne_H3_ref_esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
python ESM-2_ft/7Texas_H5_ref_esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
python ESM-2_ft/8PlosOne_H7_ref_esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
python ESM-2_ft/9PlosOne_H9_ref_esm2_mutation_effect_masking_JSY_ft1_virus10w_freeze32_lr3e-5.py
```

## vBERT-based adaptation classifier establishment, benchmarking and calculation of adaptation probability
1. H1, H3, H5, H7, and H9 sequences (segmented HA sequence consistent with the reference scheme instead of segmented HA sequences merely from 70 to 290) were sampled based on CD-HIT with a threshold of 0.98. The sampled sequences were saved as 'vBERT_classifier_establishment/H13579_PlosOne_aligned_and_sampled_sequences_without_empty_98.fasta' and 'vBERT_classifier_establishment/H13579_PlosOne_aligned_and_sampled_sequences_without_empty_98_with_origin_seqid_and_seq.csv'.

2. Hidden Markov Model-based sequence tokenization was performed on the sampled HA sequences and the vBERT embeddings of each tokenized sequence was extracted (reference: PMID: 41035817). The vBERT embeddings were saved as 'vBERT_classifier_establishment/H13579_PlosOne_aligned_and_sampled_sequences_without_empty_98_with_origin_seqid_and_seq.npy'.

3. Establishment and validation of vBERT-based adaptation classifier and quantification of the adaptation probability of each sequence in the validation dataset
```bash
python vBERT_classifier_establishment/1train_H13579_PlosOne_aligned_cdhit98_seqs_predict_prob_for_mutation_analysis_with_cm_and_roc.py
```

4. Sampling of the sequences in the validation dataset for phylogenetic analysis (MEGA was utilized for phylogenetic analysis, and the phylogenetic tree was saved as 'vBERT_classifier_establishment/phylo.nwk')
```bash
python vBERT_classifier_establishment/2create_fasta_for_phylo_sampling.py
```

5. Embedding benchmarking of vBERT and ESM-2 and ESM-2 ft on HA and Spike RBD test datasets(reference: PMID: 41035817) (vBERT and ESM-2 embeddings were evaluated on the two datasets, and we need to evaluate the fine-tuned ESM-2 on the two datasets)
```bash
python ESM-2_ft/10extract_embedding_ESM2_ft_testset_HA_JSY.py
python ESM-2_ft/10extract_embedding_ESM2_ft_testset_S_RBD_JSY.py
python ESM-2_ft/11_1padding_dimensional_reduction_clustering_HA2000_serotype_esm2_ft_coevo.py
python ESM-2_ft/11_2padding_clustering_S_RBD_1000testset_esm2_ft_coevo.py
python ESM-2_ft/12vBERT_esm2ft_esm2_embedding_benchmarking_plot.py
```
 
## CEA+ESM-2_ft / vBERT based prediction of human-adapted hidden genotypes on IAV HA sequences
1. Site importance calculation based on CEA for site 70 to 290 on HA sequences
```bash
python integrated_selection/1_0key_site_select_coevo_new_matrix_for_HA_H13579_PlosOne_aligned.py
```

2. Coupled epistatic sites calculation with RBD sites based on CEA
```bash
python integrated_selection/1_1coevo_new_matrix_for_HA_PlosOne_aligned.py
```

3. Coupled epistatic matrix calculation of site 70 to 290 on HA sequences based on CEA
```bash
python integrated_selection/1_2coevo_new_matrix_for_HA_PlosOne_aligned_with_matrix_for_plot.py
python integrated_selection/1_3coevo_new_matrix_for_HA_H135_PlosOne_aligned_plot.py
```

4. Key RBD sites selection and pairwise mutation marking with the secondary structure based on CEA
```bash
python integrated_selection/1_4site_coevo_importance_with_plot.py
python integrated_selection/1_5selected_site_and_top1_coevo_site_AApair_frequency_with_plot.py
python integrated_selection/1_6H135_selected_site1_and_coevo_top3_site2_ss_heatmap_calculation_with_plot.py
```

5. Selection top-3 coupled epistatic sites of key RBD sites and pairwise mutations analysis
```bash
python integrated_selection/2generate_mutations_with_top_sites_and_top_coevo_sites_H1.py
python integrated_selection/2generate_mutations_with_top_sites_and_top_coevo_sites_H3.py
python integrated_selection/2generate_mutations_with_top_sites_and_top_coevo_sites_H5.py
python integrated_selection/2generate_mutations_with_top_sites_and_top_coevo_sites_H7.py
python integrated_selection/2generate_mutations_with_top_sites_and_top_coevo_sites_H9.py
python integrated_selection/3H135_selected_site1_and_coevo_top3_site2_ss_heatmap_calculation_with_plot.py
```

6. Pairwise mutation re-scoring based on fine-tuned ESM-2 and vBERT-based classifier
```bash
python integrated_selection/4ESM_logist_select_H1_CA07_coevo_mutations.py
python integrated_selection/4ESM_logist_select_H3_1968_coevo_mutations.py
python integrated_selection/4ESM_logist_select_H5_Texas_coevo_mutations.py
python integrated_selection/4ESM_logist_select_H7_PlosOne_coevo_mutations.py
python integrated_selection/4ESM_logist_select_H9_PlosOne_coevo_mutations.py
python integrated_selection/5vBERT_predict_select_H1_CA07_coevo_mutations.py
python integrated_selection/5vBERT_predict_select_H3_PlosOne_coevo_mutations.py
python integrated_selection/5vBERT_predict_select_H5_Texas_coevo_mutations.py
python integrated_selection/5vBERT_predict_select_H7_PlosOne_coevo_mutations.py
python integrated_selection/5vBERT_predict_select_H9_PlosOne_coevo_mutations.py
```

7. Selection and visualization of high-mutation-effect coupled epistatic mutations based on fine-tuned ESM-2
```bash
python integrated_selection/6ESM_site12_mutation_score_freq_hist_plot.py
python integrated_selection/7ESM_select_coevo_sites_H1_plot.py
python integrated_selection/7ESM_select_coevo_sites_H3_plot.py
python integrated_selection/7ESM_select_coevo_sites_H5_plot.py
python integrated_selection/7ESM_select_coevo_sites_H7_plot.py
python integrated_selection/7ESM_select_coevo_sites_H9_plot.py
```

8. Prediction and visualization of high-adaptation coupled epistatic mutations based on vBERT classifier
```bash
python integrated_selection/8ESM_selected_before_and_after_vBERT model_prob_freq_hist_plot.py
python integrated_selection/8ESM_selected_before_and_after_vBERT model_prob_freq_hist_plot_thresh_compare_50_40_30already_100.py
python integrated_selection/9vBERT_select_coevo_sites_H1_plot.py
python integrated_selection/9vBERT_select_coevo_sites_H3_plot.py
python integrated_selection/9vBERT_select_coevo_sites_H5_plot.py
python integrated_selection/9vBERT_select_coevo_sites_H7_plot.py
python integrated_selection/9vBERT_select_coevo_sites_H9_plot.py
```

9. Mutation analysis of the final predicted site pairs (Selected pairwise mutations were saved in the 'integrated_selection/final_selected' folder)
```bash
python integrated_selection/10mutation_for_origin_software_xiantu_plot.py
python integrated_selection/11mutation_seq_distribution_and_H13_sankey_for_plot.py
python integrated_selection/11mutation_seq_distribution_and_H579_sankey_for_plot.py
python integrated_selection/12mutation_seq_distribution_and_H579_for_quantify_plot_continent.py
python integrated_selection/12mutation_seq_distribution_and_H579_for_quantify_plot_sero.py
python integrated_selection/12mutation_seq_distribution_and_H579_for_quantify_plot_sero_year.py
```

10. Literature validation of the final selected sites（The sites mentioned in the RBS-related references were saved as 'integrated_selection/literature_validation/reference_RBS_all.csv')
```bash
python integrated_selection/literature_validation/1scatterpie_references_RBS.py
python integrated_selection/literature_validation/2freq_cutoff_and_coverage_with_plot.py
```
