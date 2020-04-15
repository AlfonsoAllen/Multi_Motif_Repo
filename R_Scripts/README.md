# R SCRIPTS

* **Phenology_Week**
    * Input:
        * Raw_data/Metadata_Pollinators_Abundances_Seeds_2019.csv
    * Output:
        * Processed_data/Phenology/phenology_2019.csv: Input data + a new column that points out the week of the observed interactions
        * Processed_data/Phenology/estimated_phenology_plot.pdf: Figure that shows the estimated phenology of Caracoles by plot
        * Processed_data/Phenology/estimated_phenology_Caracoles.pdf: Figure that shows the estimated phenology of Caracoles

* **Motifs_GENERATE_DATA_WEEK_SPECIES**
    * Input:
        * Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv
    * Output:
        * Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv: Dataframe with information on the total amount of heterospecific and homospecific triplets per plot, subplot, plant species, and week for Caracoles (2019). The **triplets under consideration** in these calculations **contain two plant nodes**.

* **CORRELATIONS_motif_fitness_WEEK**
    * Input:
        * Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019.csv
        * Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv
    * Output:
        * Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/: This folder contains several figures. They show the **relation between the number of seeds and total (aggregated) amount of heterospecific and homospecific triplets** for focal plant individuals. This program **does take into account phenological co-occurrrence of plant species (by estimating triplets weekly)**.
