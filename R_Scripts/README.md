# R SCRIPTS

* **Motifs_GENERATE_DATA_WEEK**
    * Input:
        * Raw_data/Metadata_Pollinators_Abundances_Seeds_2019.csv
    * Output:
        * Processed_data/Phenology/phenology_2019.csv: Input data + a new column that points out the week of the observed interactions
        * Processed_data/Phenology/estimated_phenology_plot.pdf: Figure that shows the estimated phenology of Caracoles by plot
        * Processed_data/Phenology/estimated_phenology_Caracoles.pdf: Figure that shows the estimated phenology of Caracoles

* **Motifs_GENERATE_DATA_WEEK**
    * Input:
        *Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019.csv
    * Output:
        * Processed_data/Motifs_WEEK/Caracoles_WEEK.csv: Dataframe with information on the total amount of heterospecific and homospecific triplets per plot, subplot, plant species, and week for Caracoles (2019)

* **Motifs_GENERATE_DATA_YEAR**
    * Input:
        * Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019.csv
        * Raw_Data/abundances_2019.csv
    * Output:
        * Processed_data/Motifs_WEEK/Caracoles_WEEK.csv: Dataframe with information on the total amount of heterospecific and homospecific triplets per plot, subplot, and plant species for Caracoles (2019)
        * Processed_data/Motifs_WEEK/random_X_WEEK.csv: 100 dataframes with information on the total amount of heterospecific and homospecific triplets per plot, subplot, and plant species for Caracoles' null model.

