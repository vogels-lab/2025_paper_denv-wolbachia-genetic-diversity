############ logistic regression model ##########


Dataset_final <-  readxl::read_xlsx("DENV1-4_extractedPCR_Clean.xlsx", sheet = "Sheet1")

dissemination_rate <- Dataset_final %>% 
  mutate(`Infection(<38)` = as.numeric(`Infection(<38)`)) %>% 
  mutate(`Dissemination(<38)` = as.numeric(`Dissemination(<38)`)) %>% 
  group_by(Virus, Mosquito_strain, Serotype) %>% 
  summarise(infected = sum(`Dissemination(<38)`),
            sample_size = n())  %>% 
  ungroup() 


CT_stock <- Dataset_final %>% 
  select(Virus, `GE/uL_stock`, Serotype, Genotype) %>% 
  unique()


combined_dissemination_rate_stock <- merge(dissemination_rate, 
                                           CT_stock,
                                           by = c("Virus", "Serotype"))

combined_dissemination_rate_stock$Mosquito_strain <- factor(combined_dissemination_rate_stock$Mosquito_strain, 
                                                            levels = c("WT", "wAlbB", "wMel"))
combined_dissemination_rate_stock$Serotype <- factor(combined_dissemination_rate_stock$Serotype)
#Fit a logistic regression model
model_diss <- glm(cbind(infected, sample_size - infected) ~ log10(`GE/uL_stock`) * Mosquito_strain + Serotype,
                  family = binomial(link = "logit"),
                  data = combined_dissemination_rate_stock)
(summary(model_diss))

xtable(summary(model_diss))




# # Check for nonlinearity
# library(splines)
# model_diss_spline <- glm(
#   cbind(infected, sample_size - infected) ~ ns(log10(`GE/uL_stock`), df = 1) * Mosquito_strain + Serotype,
#   data = combined_dissemination_rate_stock,
#   family = binomial
# )
# summary(model_diss_spline)
# 
# 
# # AIC to see if the nonlinear model improves fit: lower AIC = better fit
# AIC(model_inf, model_diss_spline)



#Plot the fitted relationship
library(ggplot2)
mosquito_strain_colors <- c("WT" = "#33b89a",  
                            "wMel" = "#3b5323", 
                            "wAlbB" = "#b8cc00"
)

combined_dissemination_rate_stock$pred <- predict(model_diss, type = "response")
combined_dissemination_rate_stock$Serotype <- paste0("DENV-", combined_dissemination_rate_stock$Serotype)

Fig4 <- ggplot(combined_dissemination_rate_stock,
               aes(x = log10(`GE/uL_stock`),
                   y = 100 * infected / sample_size)) +
  geom_point(aes(fill = Mosquito_strain), color = "black", size = 3.5, shape = 21) +
  geom_line(aes(y = pred * 100, color = Mosquito_strain), size = 1.4) +
  facet_wrap(~ Serotype, nrow = 1) +
  scale_fill_manual(values = mosquito_strain_colors, limits = c("WT", "wAlbB", "wMel"), labels = c("WT", "wAlbB", "wMelM")) +  # Customize fill colors for each mosquito strain
  scale_color_manual(values = mosquito_strain_colors, limits = c("WT", "wAlbB", "wMel"), labels = c("WT", "wAlbB", "wMelM")) + # Customize colors for the smoothed lines
  theme_classic() +
  ylab("Dissemination rate (%)") + 
  xlab("Input stock virus GE/uL (log10)") + 
  guides(fill = guide_legend(title = "Mosquito colony"), color = "none") +
  theme(
    legend.position = "top",
    legend.text = element_text(color = "black", size = 20, face="bold"),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    strip.text = element_text(color = "black", size = 20, face="bold"),
    axis.text.x = element_text(color = "black", size = 20, face="bold"),
    axis.text.y = element_text(color = "black", size = 20, face="bold"),
    axis.title.x = element_text(color = "black", size = 20, face="bold"),
    axis.title.y = element_text(color = "black", size = 20, face="bold"),
    axis.line.x = element_line(color="black", linewidth = 1),
    axis.line.y = element_line(color="black", linewidth = 1),
    plot.title = element_text(size = 20, face="bold", hjust = 0.5)
  ) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) 

print(Fig4)
