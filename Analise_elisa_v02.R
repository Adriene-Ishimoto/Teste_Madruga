#install.packages("ggsignif")

library(readxl)
library(ggplot2)
library(dplyr)
library(ggsignif)

setwd("/home/Ishimoto/Documentos/PhD/DebianPC/SARS-CoV-2_qRT-PCR_Variantes/Experimentos/Elisa")

soroteca <- read_excel("/home/Ishimoto/Documentos/PhD/Soroteca_Pacientes_COVID19/Dados Brutos Soroteca/SOROTECA COVID-19_02 DE JULHO 2020_V24.xlsx")
soroteca_txt <- read.table("/home/Ishimoto/Documentos/PhD/Soroteca_Pacientes_COVID19/Dados Brutos Soroteca/SOROTECA COVID-19_02 DE JULHO 2020_V24.txt", sep = "\t")
cba_pac_tams <- read_excel("/home/Ishimoto/Documentos/PhD/Soroteca_Pacientes_COVID19/Dados Brutos CBA/CBA Total.xlsx")
cba_pac_tams$Casp1 <- as.numeric(soroteca$`ELISA Casp1 p20 (pg/ml)`[which(soroteca$Soroteca %in% cba_pac_tams$Samples)])
cba_pac_tams$IL18 <- as.numeric(soroteca$`ELISA IL-18 (pg/ml)`[which(soroteca$Soroteca %in% cba_pac_tams$Samples)])

pac_adri <- read_excel("/home/Ishimoto/Documentos/PhD/DebianPC/SARS-CoV-2_qRT-PCR_Variantes/Experimentos/Elisa/Pacientes_COVID19_Adri.xlsx")

pac_tams <- data.frame(DZ = NA, RegistroHC = soroteca$REGISTRO[which(soroteca$Soroteca %in% cba_pac_tams$Samples)],
                       Soroteca = cba_pac_tams$Samples, Nome = NA, Cepa = "COVID19_2020",
                       Casp1 = cba_pac_tams$Casp1, IL18 = cba_pac_tams$IL18, IFNy = cba_pac_tams$IFNy,
                       IL2 = cba_pac_tams$`IL-2`, IL4 = cba_pac_tams$`IL-4`, IL6 = cba_pac_tams$`IL-6`,
                       IL10 = cba_pac_tams$`IL-10`, IL17A = cba_pac_tams$`IL-17A`, TNF = cba_pac_tams$TNF)

colnames(pac_adri) <- colnames(pac_tams)

pac_adri$Casp1_adjst <- pac_adri$Casp1*3
pac_adri$IL18_adjst <- pac_adri$IL18*3
pac_adri$IL2[which(pac_adri$IL2 == 0)] <- NA
pac_adri$IL4[which(pac_adri$IL4 == 0)] <- NA
pac_adri$IL6[which(pac_adri$IL6 == 0)] <- NA
pac_adri$IL10[which(pac_adri$IL10 == 0)] <- NA
pac_adri$IL17A[which(pac_adri$IL17A == 0)] <- NA
pac_adri$IFNy[which(pac_adri$IFNy == 0)] <- NA
pac_adri$TNF[which(pac_adri$TNF == 0)] <- NA

pac_tams$IL2[which(pac_tams$IL2 == 0)] <- NA
pac_tams$IL4[which(pac_tams$IL4 == 0)] <- NA
pac_tams$IL6[which(pac_tams$IL6 == 0)] <- NA
pac_tams$IL10[which(pac_tams$IL10 == 0)] <- NA
pac_tams$IL17A[which(pac_tams$IL17A == 0)] <- NA
pac_tams$IFNy[which(pac_tams$IFNy == 0)] <- NA
pac_tams$TNF[which(pac_tams$TNF == 0)] <- NA

pac_adri$Casp1_Log <- log10(pac_adri$Casp1_adjst)
pac_adri$IL18_Log <- log10(pac_adri$IL18_adjst)
pac_adri$IFNy_Log <- log10(pac_adri$IFNy)
pac_adri$IL2_Log <- log10(pac_adri$IL2)
pac_adri$IL4_Log <- log10(pac_adri$IL4)
pac_adri$IL6_Log <- log10(pac_adri$IL6)
pac_adri$IL10_Log <- log10(pac_adri$IL10)
pac_adri$IL17A_Log <- log10(pac_adri$IL17A)
pac_adri$TNF_Log <- log10(pac_adri$TNF)

pac_tams$Casp1_Log <- log10(pac_tams$Casp1)
pac_tams$IL18_Log <- log10(pac_tams$IL18)
pac_tams$IFNy_Log <- log10(pac_tams$IFNy)
pac_tams$IL2_Log <- log10(pac_tams$IL2)
pac_tams$IL4_Log <- log10(pac_tams$IL4)
pac_tams$IL6_Log <- log10(pac_tams$IL6)
pac_tams$IL10_Log <- log10(pac_tams$IL10)
pac_tams$IL17A_Log <- log10(pac_tams$IL17A)
pac_tams$TNF_Log <- log10(pac_tams$TNF)

shapiro.test(pac_adri$Casp1_adjst)
shapiro.test(pac_adri$IL18_adjst)
shapiro.test(pac_adri$IFNy)
shapiro.test(pac_adri$IL2)
shapiro.test(pac_adri$IL4)
shapiro.test(pac_adri$IL6)
shapiro.test(pac_adri$IL10) ## normal
shapiro.test(pac_adri$IL17A)
shapiro.test(pac_adri$TNF)

shapiro.test(pac_adri$Casp1_Log)
shapiro.test(pac_adri$IL18_Log) ## normal
shapiro.test(pac_adri$IL2_Log)
shapiro.test(pac_adri$IL4_Log)
shapiro.test(pac_adri$IL6_Log)
shapiro.test(pac_adri$IL10_Log) ## normal
shapiro.test(pac_adri$IL17A_Log)
shapiro.test(pac_adri$IFNy_Log)
shapiro.test(pac_adri$TNF_Log) ## normal

shapiro.test(pac_tams$Casp1)
shapiro.test(pac_tams$IL18)
shapiro.test(pac_tams$IFNy)
shapiro.test(pac_tams$IL2) ## normal
shapiro.test(pac_tams$IL4)
shapiro.test(pac_tams$IL6)
shapiro.test(pac_tams$IL10)
shapiro.test(pac_tams$IL17A)
shapiro.test(pac_tams$TNF)

shapiro.test(pac_tams$Casp1_Log)
shapiro.test(pac_tams$IL18_Log) ## normal
shapiro.test(pac_tams$IL2_Log) ## normal
shapiro.test(pac_tams$IL4_Log) ## normal
shapiro.test(pac_tams$IL6_Log) ## normal
shapiro.test(pac_tams$IL10_Log) ## normal
shapiro.test(pac_tams$IL17A_Log)
shapiro.test(pac_tams$IFNy_Log)
shapiro.test(pac_tams$TNF_Log) ## normal

pac_adri_tams <- rbind(pac_adri[, -c(15,16)], pac_tams)

pac.casp1.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(Casp1_Log, na.rm = TRUE),
    Casp1_Log = mean(Casp1_Log, na.rm = TRUE)
  )

pac.il18.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL18_Log, na.rm = TRUE),
    IL18_Log = mean(IL18_Log, na.rm = TRUE)
  )

pac.il2.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL2_Log, na.rm = TRUE),
    IL2_Log = mean(IL2_Log, na.rm = TRUE)
  )

pac.il4.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL4_Log, na.rm = TRUE),
    IL4_Log = mean(IL4_Log, na.rm = TRUE)
  )

pac.il6.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL6_Log, na.rm = TRUE),
    IL6_Log = mean(IL6_Log, na.rm = TRUE)
  )

pac.il10.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL10_Log, na.rm = TRUE),
    IL10_Log = mean(IL10_Log, na.rm = TRUE)
  )

pac.il17a.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IL17A_Log, na.rm = TRUE),
    IL17A_Log = mean(IL17A_Log, na.rm = TRUE)
  )

pac.ifny.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(IFNy_Log, na.rm = TRUE),
    IFNy_Log = mean(IFNy_Log, na.rm = TRUE)
  )

pac.tnf.summary <- pac_adri_tams %>%
  group_by(Cepa) %>%
  summarise(
    sd = sd(TNF_Log, na.rm = TRUE),
    TNF_Log = mean(TNF_Log, na.rm = TRUE)
  )


casp1 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = Casp1_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.casp1.summary[-which(is.na(pac.casp1.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = Casp1_Log-sd, ymax = Casp1_Log+sd), data = pac.casp1.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("Casp-1 p20") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL18 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL18_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il18.summary[-which(is.na(pac.il18.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL18_Log-sd, ymax = IL18_Log+sd), data = pac.il18.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-18") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL2 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL2_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il2.summary[-which(is.na(pac.il2.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL2_Log-sd, ymax = IL2_Log+sd), data = pac.il2.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-2") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL4 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL4_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il4.summary[-which(is.na(pac.il4.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL4_Log-sd, ymax = IL4_Log+sd), data = pac.il4.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-4") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL6 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL6_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il6.summary[-which(is.na(pac.il6.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL6_Log-sd, ymax = IL6_Log+sd), data = pac.il6.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-6") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL10 <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL10_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il10.summary[-which(is.na(pac.il10.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL10_Log-sd, ymax = IL10_Log+sd), data = pac.il10.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-10") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IL17A <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IL17A_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.il17a.summary[-which(is.na(pac.il17a.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IL17A_Log-sd, ymax = IL17A_Log+sd), data = pac.il17a.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IL-17A") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

IFNy <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = IFNy_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.ifny.summary[-which(is.na(pac.ifny.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = IFNy_Log-sd, ymax = IFNy_Log+sd), data = pac.ifny.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("IFNy") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

TNF <- ggplot(pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(x = Cepa, y = TNF_Log, fill = Cepa)) +
  geom_bar(stat = "identity", data = pac.tnf.summary[-which(is.na(pac.tnf.summary$Cepa)),], 
           color = "black", width = 0.6, show.legend = FALSE) +
  geom_jitter(data = pac_adri_tams[-which(is.na(pac_adri_tams$Cepa)),], aes(shape = Cepa), position = position_jitter(0.2), 
              color = "black", size = 4, show.legend = FALSE) + 
  geom_errorbar(aes(ymin = TNF_Log-sd, ymax = TNF_Log+sd), data = pac.tnf.summary, width = 0.3) +
  scale_fill_manual(values = c("lightcoral", "lightblue1", "dodgerblue3", "steelblue1")) +
  scale_x_discrete(limits = c("COVID19_2020", "CT", "WT", "P1"), labels = c("COVID19_2020" = "COVID-19 (2020)", "CT" = "Healthy Donors", "WT" = "Parental", "P1" = "Gamma")) + 
  geom_signif(comparisons = list(c("CT", "WT"), c("CT", "P1"), c("WT", "P1"), c("CT", "COVID19_2020"), 
                                 c("COVID19_2020", "WT"), c("COVID19_2020", "P1")), map_signif_level = TRUE,
              y_position = c(4.8, 4.4, 4.1, 4, 5.2, 5.5), size = 0.6, textsize = 7) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) + 
  ggtitle("TNF") + 
  xlab("") + 
  ylab("Log10 (pg/mL)") + 
  theme(plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
        axis.text = element_text(size = 18, colour = "black"), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 

#############################################################################################

ggsave("Casp1_v02.png", casp1, dpi = 300, height = 9, width = 9)
ggsave("IL18_v02.png", IL18, dpi = 300, height = 9, width = 9)
ggsave("IL2.png", IL2, dpi = 300, height = 9, width = 9)
ggsave("IL4.png", IL4, dpi = 300, height = 9, width = 9)
ggsave("IL6.png", IL6, dpi = 300, height = 9, width = 9)
ggsave("IL10.png", IL10, dpi = 300, height = 9, width = 9)
ggsave("IL17A.png", IL17A, dpi = 300, height = 9, width = 9)
ggsave("IFNy.png", IFNy, dpi = 300, height = 9, width = 9)
ggsave("TNF.png", TNF, dpi = 300, height = 9, width = 9)