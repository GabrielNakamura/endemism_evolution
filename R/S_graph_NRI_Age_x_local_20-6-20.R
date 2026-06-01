# Performing unpaired (two independent groups) analysis.
unpaired_mean_diff_age<- dabestr::dabest(anova_data, local, mean_age,
                             idx = c("temperate", "tropical"),
                             paired = FALSE)
unpaired_mean_diff_NRI<- dabestr::dabest(anova_data, local, NRI,
                                         idx = c("temperate", "tropical"),
                                         paired = FALSE)
unpaired_mean_diff_sesMpd<- dabestr::dabest(anova_data_ses, local, NRI,
                                         idx = c("temperate", "tropical"),
                                         paired = FALSE)

# Display the results in a user-friendly format.
quartz()
plot(unpaired_mean_diff_age)
plot(unpaired_mean_diff_NRI)
plot(unpaired_mean_diff_sesMpd)


