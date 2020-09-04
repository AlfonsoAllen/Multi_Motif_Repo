# Dot chart of a single numeric vector

fitness.data_ho <- fitness.data[order(fitness.data$homo_motif), ]
dotchart(fitness.data_ho$homo_motif, labels = row.names(paste0(
  fitness.data_ho$Plot," ",fitness.data_ho$Subplot," ",fitness.data_ho$Plant_Simple)),
  cex = 0.6, xlab = "homo_motifs")

fitness.data_he <- fitness.data[order(fitness.data$hete_motif), ]
dotchart(fitness.data_he$hete_motif, labels = row.names(paste0(
  fitness.data_he$Plot," ",fitness.data_he$Subplot," ",fitness.data_he$Plant_Simple)),
  cex = 0.6, xlab = "hete_motif")

fitness.data_st <- fitness.data[order(fitness.data$StrengthIn), ]
dotchart(fitness.data_st$StrengthIn, labels = row.names(paste0(
  fitness.data_st$Plot," ",fitness.data_st$Subplot," ",fitness.data_st$Plant_Simple)),
  cex = 0.6, xlab = "StrengthIn")

fitness.data_r <- fitness.data[order(fitness.data$Ratio), ]
dotchart(fitness.data_r$Ratio, labels = row.names(paste0(
  fitness.data_r$Plot," ",fitness.data_r$Subplot," ",fitness.data_r$Plant_Simple)),
  cex = 0.6, xlab = "Ratio")

