
############Convert data to relative abundances#################

norm.data <- read.table("Combined_Study4_DeSeq2_normalised_data_table.txt")
dim(norm.data)
norm.data <- as.matrix(norm.data)

metadata <- read.table("Combined_Study4_metadata.txt", sep = "\t", header = TRUE, row.names = 1)
summary(metadata$FLG_mutation)

rel.data <- (make_relative(norm.data)*100)
dim(rel.data)
#rowSums(rel.data)
#data <- data[, colSums(data != 0) > 0]

row.names(metadata) <- row.names(rel.data)
attach(metadata)
dim(metadata)


Taxonomy <- read.table("Combined_Stud1_2_and_4_taxonomy.txt", sep = '\t', header = TRUE)
row.names(Taxonomy) <- colnames(rel.data)
colnames(Taxonomy) <- "taxa"
Taxonomy <- str_split_fixed(Taxonomy$taxa, ";", 6)
rownames(Taxonomy) <- colnames(norm.data)
colnames(Taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
dim(Taxonomy)

#####Remove ASVs that are in high occurrence wihtin the negatives
#####Less than 10.

NegData <- rel.data[c(which(metadata$Individual == "Negative")),]
dim(NegData)
NegMetadata <- metadata[c(which(metadata$Individual == "Negative")),]
dim(NegMetadata)
NegASVPersistences <- colSums(NegData > 0)
rel.data <- rel.data[, -c(which(NegASVPersistences > 10))]
dim(rel.data)
mean(rowSums(rel.data))
Taxonomy <- Taxonomy[-c(which(NegASVPersistences > 10)), ]
dim(Taxonomy)
Taxonomy <- as.data.frame(Taxonomy)

rel.data <- rel.data[-c(which(metadata$Individual == "Negative")),]
dim(rel.data)
metadata <- metadata[-c(which(metadata$Individual == "Negative")),]
dim(metadata)
#Need to find a taxonomy file to go with this. 


Archaea <-  which(Taxonomy$Kingdom == "Archaea")
mean(rowSums(rel.data[,c(Archaea)]))
Eukaryota <- which(Taxonomy$Kingdom == "Eukaryota")
mean(rowSums(rel.data[,c(Eukaryota)]))
NA.taxa <- which(Taxonomy$Kingdom == "NA")
mean(rowSums(rel.data[, NA.taxa]))
Bacteria.taxa <- which(Taxonomy$Kingdom == "Bacteria")
mean(rowSums(rel.data[, Bacteria.taxa]))

Not.bacteria.taxa <- which(Taxonomy$Kingdom != "Bacteria")
mean(rowSums(rel.data[, Not.bacteria.taxa]))

rel.data <- rel.data[, -c(Not.bacteria.taxa)]
Taxonomy <- Taxonomy[-c(Not.bacteria.taxa), ]
dim(Taxonomy)

low.read.samples <- which(rowSums(rel.data) < 50)

rel.data <- rel.data[-c(low.read.samples),]
dim(rel.data)
metadata <- metadata[-c(low.read.samples),]
dim(metadata)

Tape10Samples <- which(metadata$TapeNo == 10)
rel.data <- rel.data[-c(Tape10Samples),]
dim(rel.data)
metadata <- metadata[-c(Tape10Samples),]
dim(metadata)

rel.data <- (make_relative(rel.data)*100)
dim(rel.data)

Richness <- rowSums(rel.data > 0)
metadata <- cbind(metadata, Richness)

ShannonDiv <- diversity(as.data.frame(rel.data), index = "shannon")
metadata <- cbind(metadata, ShannonDiv)

metadata$Visit <- as.factor(metadata$Visit)
metadata$Individual <- as.factor(metadata$Individual)
metadata$VisitLinear <- as.integer(droplevels(metadata$Visit))


#rowSums(rel.data)
#data <- data[, colSums(data != 0) > 0]

row.names(metadata) <- row.names(rel.data)

metadata$PatientID <- factor(metadata$PatientID, levels = c(levels =
  "AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD7", "AD8", "AD9", "AD10", 
  "AD11", "AD12", "AD13", "AD14", "AD15", "AD16", "AD17", "AD18", "AD19", "AD20", 
  "HC1", "HC2", "HC3", "HC4", "HC5", "HC6", "HC7", "HC8", "HC9", "HC10", 
  "HC11", "HC12", "HC13", "HC14", "HC15", "HC16", "HC17", "HC18", "HC19", "HC20",
  "HC21", "HC22", "HC23"),
  labels = c(
  "AD1", "AD2", "AD3", "AD4", "AD5", "AD6", "AD7", "AD8", "AD9", "AD10", 
  "AD11", "AD12", "AD13", "AD14", "AD15", "AD16", "AD17", "AD18", "AD19", "AD20", 
  "HC1", "HC2", "HC3", "HC4", "HC5", "HC6", "HC7", "HC8", "HC9", "HC10", 
  "HC11", "HC12", "HC13", "HC14", "HC15", "HC16", "HC17", "HC18", "HC19", "HC20",
  "HC21", "HC22", "HC23")
)

metadata$TapeNo <- ordered(metadata$TapeNo, levels = c("1", "5", "10", "Negative"),
                           labels = c("Tape 1", "Tape 5", "Tape 10", "Negative"))

metadata$pH <- gsub("missing", NA, metadata$pH)
metadata$pH <- as.numeric(as.character(metadata$pH))

metadata$O.SCORAD <- gsub("<NA>", NA, metadata$O.SCORAD)
metadata$O.SCORAD <- as.numeric(as.character(metadata$O.SCORAD))

metadata$TEWL <- gsub("missing", NA, metadata$TEWL)
metadata$TEWL <- as.numeric(as.character(metadata$TEWL))

metadata$SteroidUse <- factor(metadata$SteroidUse, levels = c(levels =
  "Negative", "None", "Short", "Mid", "Long"),
  labels = c("Negative", "None", "Short", "Mid", "Long")
)


brewer.pal(n = 3, name = "Dark2")
TEWL.plot <- 
  ggplot(metadata, aes(x = PatientID, y = TEWL, 
                       fill = Individual)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_grid( ~ Individual, scales = "free_x") +
  facet_grid( ~ Individual, scales = "free_x") + labs(x = NULL, fill = NULL) + 
  scale_fill_manual(values = c(brewer.pal(n = 3, name = "Set2"))) + ggtitle("b") +
  theme(axis.text.x = element_text(angle = 90, size = 8))
TEWL.plot

ggplot(metadata, aes(x = (sqrt(TEWL)), y = (sqrt(Richness)))) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  stat_smooth(data = metadata, method = "lm", formula = y ~ x, se = FALSE) + facet_grid( ~ Individual, scales = "free_x") 

pH.plot <- ggplot(metadata, aes(x = PatientID, y = pH, 
                                fill = Individual)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_grid( ~ Individual, scales = "free_x") +
  ggtitle("a") + labs(x = NULL, fill = NULL) + 
  scale_fill_manual(values = c(brewer.pal(n = 3, name = "Set2"))) +
  theme(axis.text.x = element_text(angle = 90, size = 8))
pH.plot

ggplot(metadata, aes(x = pH, y = Richness)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  stat_smooth(data = metadata, method = "lm", formula = y ~ x, se = FALSE) + facet_grid( ~ Individual, scales = "free_x") 

ggplot(metadata, aes(x = O.SCORAD, y = Richness)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  stat_smooth(data = metadata, method = "lm", formula = y ~ x, se = FALSE) #+ facet_grid( ~ Individual, scales = "free_x") 

SCORAD.plot <- 
  ggplot(metadata, aes(x = PatientID, y = O.SCORAD, 
  fill = Individual)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_wrap( ~ Individual, scales = "free_x", drop = TRUE) +
  ggtitle("c") + labs(x = NULL, fill = NULL, y = "SCORAD") + 
  scale_fill_manual(values = c(brewer.pal(n = 3, name = "Set2"))) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 8))
SCORAD.plot

ggplot(metadata, aes(x = O.SCORAD, y = Richness)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
  stat_smooth(data = metadata, method = "lm", formula = y ~ x, se = FALSE) #+ facet_grid( ~ Individual, scales = "free_x") 


ggplot(metadata, aes(x = PatientID, y = Richness, 
                     fill = SteroidUse)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_grid( ~ Individual, scales = "free_x") 

ggplot(metadata, aes(x = PatientID, y = Richness, 
                     fill = TopicalAB)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_grid( ~ Individual, scales = "free_x") 

ggplot(metadata, aes(x = PatientID, y = Richness, 
                     fill = AB.time)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + facet_grid( ~ Individual, scales = "free_x") 


library(patchwork)
combined <- pH.plot + TEWL.plot + SCORAD.plot & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect") + plot_layout(ncol = 3)


##############ALL samples compared (without negatives)###############

x <- (colSums(rel.data != 0)) > 0 #OTU as columns
rel.data <- rel.data[, x == "TRUE"] #Removes taxa that are all 0s
dim(rel.data)

Taxonomy <- Taxonomy[ x == "TRUE",] #Removes taxa that are all 0s
dim(Taxonomy)



AllDataBray.ordination <- rel.data[-c(which(row.names(rel.data) == "Study_4.s404")),]
dim(AllDataBray.ordination)
AllMetaataBray.ordination <- metadata[-c(which(row.names(metadata) == "Study_4.s404")),]
dim(AllMetaataBray.ordination)

#AllDataBray <- vegdist(AllDataBray.ordination, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#AllDataMDS <- monoMDS(AllDataBray) # Performs multidimensional scaling for the data
plot(AllDataMDS)

Nmds1 <- AllDataMDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- AllDataMDS$points[,2] #Extracts dimension 2 as a vector

Nmds=cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds<-as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

#bioenv(AllDataBray.ordination ~ AllMetaataBray.ordination$Individual + AllMetaataBray.ordination$LS.NLS + 
#  AllMetaataBray.ordination$PatientID + AllMetaataBray.ordination$Samp_location + AllMetaataBray.ordination$TapeNo +
#    AllMetaataBray.ordination$Visit, metric = "gower", index = "bray") #Input matrix is binary

#cca(AllDataBray.ordination ~ AllMetaataBray.ordination$PatientID)


model <- glm(AllMetaataBray.ordination$Richness ~ AllMetaataBray.ordination$PatientID + 
               AllMetaataBray.ordination$LS.NLS + 
               AllMetaataBray.ordination$TapeNo + AllMetaataBray.ordination$Visit + 
               AllMetaataBray.ordination$Samp_location)
drop1(model, test = "Chisq")

model <- glm(AllMetaataBray.ordination$ShannonDiv ~ AllMetaataBray.ordination$PatientID + 
               AllMetaataBray.ordination$LS.NLS + 
               AllMetaataBray.ordination$TapeNo + 
               AllMetaataBray.ordination$Visit + AllMetaataBray.ordination$Samp_location)
drop1(model, test = "Chisq")




adonis(AllDataBray ~  AllMetaataBray.ordination$LS.NLS + 
         AllMetaataBray.ordination$Samp_location + AllMetaataBray.ordination$TapeNo +
         AllMetaataBray.ordination$Visit, strata = AllMetaataBray.ordination$PatientID)

model <- lmer(AllMetaataBray.ordination$Richness ~ AllMetaataBray.ordination$LS.NLS + 
                AllMetaataBray.ordination$Samp_location + 
                AllMetaataBray.ordination$TapeNo + AllMetaataBray.ordination$Visit + (1|AllMetaataBray.ordination$PatientID))
drop1(model, test = "Chisq")

model <- lmer(AllMetaataBray.ordination$Richness ~  AllMetaataBray.ordination$Visit + (1|AllMetaataBray.ordination$PatientID))
drop1(model, test = "Chisq")

model <- lmer(AllMetaataBray.ordination$ShannonDiv ~ AllMetaataBray.ordination$Individual + AllMetaataBray.ordination$Samp_location + 
                AllMetaataBray.ordination$TapeNo + AllMetaataBray.ordination$Visit + (1|AllMetaataBray.ordination$PatientID))
drop1(model, test = "Chisq")

model <- lmer(AllMetaataBray.ordination$ShannonDiv ~ AllMetaataBray.ordination$Visit + (1|AllMetaataBray.ordination$PatientID))
drop1(model, test = "Chisq")


sample.MDS <-ggplot(data=Nmds, aes(y=Nmds2, x=Nmds1, Type="p", colour = AllMetaataBray.ordination$PatientID, 
  shape = AllMetaataBray.ordination$LS.NLS)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2") + 
  theme(plot.title=element_text(hjust=0)) +
  #scale_shape_manual(name = "Patient ID",
  #values = c(15, 16, 17, 15, 16, 17, 15, 16, 17, 15, 16, 17,
  #15, 16, 17, 15, 16, 17, 15, 16, 17, 15, 16, 17,
  #15, 16, 17, 15, 16, 17, 15, 16, 17, 15, 16, 17,
  #15, 16, 17, 15, 16, 17, 15)) + 
  ggtitle("a") + 
  theme(text = element_text(size=12)) + 
  #theme(legend.text=element_text(size=8)) +
  scale_shape_manual(name = "Skin status", values = c(20, 18, 15)) + theme(legend.position = "none") +
  theme(text = element_text(size=12))

colourCount = length(unique(metadata$PatientID))
getPalette = colorRampPalette(brewer.pal(43, "Accent"))

Patient.MDS <- sample.MDS + scale_colour_manual(name = "Patient ID", values = getPalette(colourCount), guide = FALSE)
Patient.MDS #Shows you our graph

sample.MDS <-ggplot(data=Nmds, aes(y=Nmds2, x=Nmds1, Type="p", colour = AllMetaataBray.ordination$TapeNo,
  shape = AllMetaataBray.ordination$LS.NLS)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = "Tape number") + 
  theme(plot.title=element_text(hjust=0)) + ggtitle("b")  + theme(text = element_text(size=10)) +
  scale_shape_manual(name = "Skin status", values = c(20, 18, 15)) 
Depth.nMDS <- sample.MDS
Depth.nMDS

metadata$Label <- paste(metadata$PatientID, metadata$LS.NLS, sep="-")
metadata$Label <- str_remove(metadata$Label, "-HC")
metadata$Layer <- ordered(metadata$TapeNo, labels = c("Surface", "Within Epi."))

TapeRichnessBoxplot <- 
  ggplot(metadata, aes(x = Label, y = Richness, 
  fill = Layer)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  + 
  facet_grid( ~ Individual, scales = "free_x", space='free_x') +
  labs(x = NULL, fill = NULL) + ggtitle("a") + theme(text = element_text(size=12))
TapeRichnessBoxplot


grid.arrange(Patient.MDS, Depth.nMDS, TapeRichnessBoxplot, 
  layout_matrix = rbind(c(1, 1, 2, 2, 2), c(3,3,3,3,3)))



model <- lmer(metadata$Richness ~ metadata$TapeNo +
  (1|metadata$PatientID:metadata$LS.NLS), data = metadata)
drop1(model, test = "Chisq")

model <- lmer(metadata$ShannonDiv ~ metadata$TapeNo +
  (1|metadata$PatientID:metadata$LS.NLS), data = metadata)
drop1(model, test = "Chisq")

aggregate(metadata$Richness, by = list(Category = metadata$TapeNo), FUN = mean)
aggregate(metadata$ShannonDiv, by = list(Category = metadata$TapeNo), FUN = mean)








metadata$DepthVariable <- paste(metadata$PatientID, metadata$LS.NLS, metadata$Samp_location, metadata$Visit, sep="_")

frequency <- as.data.frame(metadata %>% group_by(DepthVariable) %>% tally())
z <- which(frequency$n == "1")
z <- frequency[c(z),1]

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
MetadataPersistences <- subset(metadata, subset = DepthVariable %not in% c(z))
dim(MetadataPersistences)
RelDataPersistences <- subset(rel.data, subset = metadata$DepthVariable %not in% c(z))
dim(RelDataPersistences)

BinaryDataPersistences <- (RelDataPersistences > 0) # logical, or
BinaryDataPersistences <- (BinaryDataPersistences > 0)*1L # integer 01

a <- which(MetadataPersistences$TapeNo == "Tape 1")
b <- which(MetadataPersistences$TapeNo == "Tape 5")

dim(MetadataPersistences)
dim(BinaryDataPersistences)

Tape1Metadata <- MetadataPersistences[a,]
Tape5Metadata <- MetadataPersistences[b,]

Tape1s <- BinaryDataPersistences[a,]
dim(Tape1s)
Tape5s <- BinaryDataPersistences[b,]
dim(Tape5s)

rownames(Tape1s) <- Tape1Metadata$DepthVariable
rownames(Tape5s) <- Tape5Metadata$DepthVariable

Tape1Metadata <- Tape1Metadata[order(Tape1Metadata$DepthVariable),]
dim(Tape1Metadata)
Tape5Metadata <- Tape5Metadata[order(Tape5Metadata$DepthVariable),]
dim(Tape5Metadata)

Tape1s <- Tape1s[order(rownames(Tape1s)),]
Tape5s <- Tape5s[order(rownames(Tape5s)),]

Persistences <- as.data.frame(Tape5s - Tape1s)
class(Persistences)
PersistenceMeans <- colMeans(Persistences)


TotalTapePersistences <- (colSums(Tape5s) + colSums(Tape1s))
TotalTapePersistenceDifferences <- (colSums(Tape5s) - colSums(Tape1s))
Tape1sReduced <- Tape1s[, c(which(TotalTapePersistences > 26))]
dim(Tape1sReduced)
Tape5sReduced <- Tape5s[, c(which(TotalTapePersistences > 26))]
dim(Tape5sReduced)
TaxonomyReduced <- Taxonomy[c(which(TotalTapePersistences > 26)), ]
dim(TaxonomyReduced)

TapesData <- rbind(Tape1sReduced, Tape5sReduced)
dim(TapesData)
TapesMetadata <- rbind(Tape1Metadata, Tape5Metadata)
dim(TapesMetadata)
TapesMetadata$TapeNo <- droplevels(TapesMetadata$TapeNo)


Tape1PatientPersistences <- aggregate(Tape1sReduced, 
  by = list(Category = Tape1Metadata$LS.NLS, Tape1Metadata$PatientID), FUN = mean)
dim(Tape1PatientPersistences)
Tape1PatientPersistencesMetadata <- Tape1PatientPersistences[ , 1:2 ]
dim(Tape1PatientPersistencesMetadata)
Tape1PatientPersistences <- Tape1PatientPersistences[ , -c(1:2) ]
dim(Tape1PatientPersistences)

Tape5PatientPersistences <- aggregate(Tape5sReduced, 
  by = list(Category = Tape5Metadata$LS.NLS, Tape5Metadata$PatientID), FUN = mean)
dim(Tape5PatientPersistences)
Tape5PatientPersistencesMetadata <- Tape5PatientPersistences[ , 1:2 ]
dim(Tape5PatientPersistencesMetadata)
Tape5PatientPersistences <- Tape5PatientPersistences[ , -c(1:2) ]
dim(Tape5PatientPersistences)

t.test.Outputs.t <- list()
t.test.Outputs.p <- list()

for (i in 1:ncol(Tape1PatientPersistences)){
  t.test.Output <- t.test(Tape1PatientPersistences[,i], Tape5PatientPersistences[,i], paired = TRUE)
  t.test.Outputs.t[[i]] <- cbind(t.test.Output[[1]])
  t.test.Outputs.p[[i]] <- cbind(t.test.Output[[3]])
}

t.test.Outputs.t <- unlist(t.test.Outputs.t, use.names = FALSE)
t.test.Outputs.p <- unlist(t.test.Outputs.p, use.names = FALSE)

t.test.Outputs <- data.frame("t-value" = t.test.Outputs.t, "P-value" = t.test.Outputs.p)
dim(t.test.Outputs)

dim(Tape1PatientPersistences)
dim(Tape5PatientPersistences)
t.test.Outputs$q.value <- p.adjust(t.test.Outputs$P.value, method = "fdr")
t.test.Outputs <- data.frame(TaxonomyReduced, t.test.Outputs)
write.table(t.test.Outputs, file = "t_tests_depth_results.txt", sep = "\t")

plot(rowSums(Tape1PatientPersistences), rowSums(Tape5PatientPersistences), xlab = "Tape 1 persistences", ylab = "Tape 5 persistences") 
abline(0, 1, lwd = 1, lty = 5)
abline(0,0.53, lwd = 1, col = "green", lty = 5)

TapesReduced <- data.frame("Tape1s" = (rowSums(Tape1PatientPersistences)), "Tape5s" = (rowSums(Tape5PatientPersistences)))
TapesReducedSignificant <- (TapesReduced[c(a), ])
ggplot(TapesReduced, aes(x = Tape1s, y = Tape5s)) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_abline(slope=1, colour = "Green") + geom_abline(slope = 0.53, linetype = 2, colour = "Red") +
  geom_point(data=TapesReducedSignificant, aes(x = Tape1s, y = Tape5s), color='red',size=3)

a <- which(t.test.Outputs$q.value < 0.05)
a <- which(t.test.Outputs$q.value < 0.00001) 
Tape1sReducedPersistancesSignificant <- Tape1PatientPersistences[, c(a)]
dim(Tape1sReducedPersistancesSignificant)
Tape5sReducedPersistancesSignificant <- Tape5PatientPersistences[, c(a)]
dim(Tape5sReducedPersistancesSignificant)
TaxonomyReducedSignificant <- TaxonomyReduced[c(a), ]
dim(TaxonomyReducedSignificant)

TapesReduced <- data.frame("Tape1s" = (colMeans(Tape1sReducedPersistancesSignificant)), 
  "Tape5s" = (colMeans(Tape5sReducedPersistancesSignificant)))
dim(TapesReduced)
ggplot(TapesReduced, aes(x = Tape1s, y = Tape5s, label = row.names(TaxonomyReducedSignificant))) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_abline(slope=1, linetype = 2, colour = "Blue") +
  geom_abline(slope = 0.53, linetype = 2, colour = "Red") +
  geom_point(data=TapesReduced, aes(x = Tape1s, y = Tape5s), color='black',size=3) + 
  labs(x = "Surface", y = "Within Epi.") #+ geom_label() + geom_label_repel()


dim(Tape1sReducedPersistancesSignificant)
Tape1PatientPersistences <- aggregate(Tape1sReducedPersistancesSignificant, 
  by = list(Category = Tape1PatientPersistencesMetadata$Category), FUN = mean)

Tape5PatientPersistences <- aggregate(Tape5sReducedPersistancesSignificant, 
  by = list(Category = Tape5PatientPersistencesMetadata$Category), FUN = mean)

Tape1PatientPersistences$TapeNo <- c(rep("Tape1", times = nrow(Tape1PatientPersistences)))
Tape5PatientPersistences$TapeNo <- c(rep("Tape5", times = nrow(Tape5PatientPersistences)))

wide.data <- rbind(Tape1PatientPersistences, Tape5PatientPersistences)
long.data <- melt(wide.data, id.var = c("Category", "TapeNo"))

Tape1PatientPersistencesSDs <- aggregate(Tape1sReducedPersistancesSignificant, 
  by = list(Category = Tape1PatientPersistencesMetadata$Category), FUN = sd)

Tape5PatientPersistencesSDs <- aggregate(Tape5sReducedPersistancesSignificant, 
  by = list(Category = Tape5PatientPersistencesMetadata$Category), FUN = sd)

Tape1PatientPersistences$TapeNo <- c(rep("Tape1", times = nrow(Tape1PatientPersistences)))
Tape5PatientPersistences$TapeNo <- c(rep("Tape5", times = nrow(Tape5PatientPersistences)))

wide.data <- rbind(Tape1PatientPersistences, Tape5PatientPersistences)
long.dataSDs <- melt(wide.data, id.var = c("Category", "TapeNo"))

long.data$SDs <- as.numeric(as.character(long.data$value))
colnames(long.data) <- c("Status", "TapeNo", "SeqVar", "Abundance", "SD")
long.data$SEM <- long.data$SD/(sqrt(nrow(long.data)))


label <- paste(TaxonomyReducedSignificant[,6], "sp.", sep=" ")
label2 <- gsub('SeqVar', 'ASV', row.names(TaxonomyReducedSignificant))
label3 <- paste("(", label2, ")", sep="")
label4 <- paste(label, label3, sep=" ")
label5 <- label4[-9]

toBeRemoved <- which(long.data$SeqVar=="SeqVar44")
long.data1 <- long.data[-toBeRemoved,]
dim(long.data1)

long.data1$Layer <- ordered(long.data1$TapeNo, labels = c("Surface", "Within Epi."))

TapePersistencesPlot <- 
  ggplot(long.data1, aes(x = SeqVar, y = Abundance, fill = Layer)) +
  geom_bar(stat = "identity", colour = "black", position = "dodge") + 
  theme(strip.text.x = element_text(angle = 0)) + theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  facet_wrap( ~ Status, scales = "fixed", ncol = 1, strip.position="right") +
  geom_errorbar(aes(ymin = Abundance - SEM, 
  ymax = Abundance + SEM), width=.2,
  position=position_dodge(.9)) + scale_x_discrete(labels=c(label5)) +
  labs(x = NULL, y = "Persistence", fill = NULL) + scale_fill_manual(values = c("#490052", "#ffe741")) +
  ggtitle("b") + guides(fill = FALSE)
TapePersistencesPlot



library(patchwork)
combined <- TapeRichnessBoxplot + TapePersistencesPlot & theme(legend.position = "right")
combined + plot_layout(guides = "collect") + plot_layout(ncol = 1)







###########Turnover data setup#############

BigDataJaccards <- readRDS(file = "BigDataTape1JaccardsData.RDS")
dim(BigDataJaccards)
BigDataJaccardsFiltered <- BigDataJaccards[-c(low.read.samples), -c(low.read.samples)]
BigDataJaccardsFiltered <- BigDataJaccardsFiltered[-c(Tape10Samples), -c(Tape10Samples)]
dim(BigDataJaccardsFiltered)

BigDataJaccardsVectors1 <- dist_groups(BigDataJaccardsFiltered, metadata$LS.NLS)
BigDataJaccardsVectors2 <- dist_groups(BigDataJaccardsFiltered, metadata$PatientID)
BigDataJaccardsVectors3 <- dist_groups(BigDataJaccardsFiltered, metadata$TapeNo)

BigdataYC <- readRDS(file = "BigdataTape1YCMatrix.RDS")
BigdataYCFiltered <- BigdataYC[-c(low.read.samples), -c(low.read.samples)]
BigdataYCFiltered <- BigdataYCFiltered[-c(Tape10Samples), -c(Tape10Samples)]
dim(BigdataYCFiltered)

BigdataYCVectors1 <- dist_groups(BigdataYCFiltered, metadata$LS.NLS)
BigdataYCVectors2 <- dist_groups(BigdataYCFiltered, metadata$PatientID)
BigdataYCVectors3 <- dist_groups(BigdataYCFiltered, metadata$TapeNo)


TimeDifferences <- outer(metadata$VisitLinear, metadata$VisitLinear, '-')
TimeDifferences <- abs(TimeDifferences)
dim(TimeDifferences)

TimeDifferencesVectors <- dist_groups(TimeDifferences, metadata$PatientID)
TimeDifferencesVectors$Distance <- abs(TimeDifferencesVectors$Distance)
TimeDifferencesVectors$Distance


pHDifferences <- outer(metadata$pH, metadata$pH, '-')
#pHDifferences <- abs(pHDifferences)
dim(pHDifferences)

pHDifferencesVectors <- dist_groups(pHDifferences, metadata$PatientID)
#pHDifferencesVectors$Distance <- abs(pHDifferencesVectors$Distance)
pHDifferencesVectors$Distance


SCORADDifferences <- outer(metadata$O.SCORAD, metadata$O.SCORAD, '-')
#SCORADDifferences <- abs(SCORADDifferences)
dim(SCORADDifferences)

SCORADDifferencesVectors <- dist_groups(SCORADDifferences, metadata$PatientID)
#SCORADDifferencesVectors$Distance <- abs(SCORADDifferencesVectors$Distance)
SCORADDifferencesVectors$Distance

TEWLDifferences <- outer(metadata$TEWL, metadata$TEWL, '-')
#TEWLDifferences <- abs(TEWLDifferences)
dim(TEWLDifferences)

TEWLDifferencesVectors <- dist_groups(TEWLDifferences, metadata$PatientID)
#TEWLDifferencesVectors$Distance <- abs(TEWLDifferencesVectors$Distance)
TEWLDifferencesVectors$Distance

RichnessDifferences <- outer(metadata$Richness, metadata$Richness, '-')
#RichnessDifferences <- abs(RichnessDifferences)
dim(RichnessDifferences)

RichnessDifferencesVectors <- dist_groups(RichnessDifferences, metadata$PatientID)
#RichnessDifferencesVectors$Distance <- abs(RichnessDifferencesVectors$Distance)
RichnessDifferencesVectors$Distance

ShannonDifferences <- outer(metadata$ShannonDiv, metadata$ShannonDiv, '-')
#ShannonDifferences <- abs(ShannonDifferences)
dim(ShannonDifferences)

ShannonDifferencesVectors <- dist_groups(ShannonDifferences, metadata$PatientID)
#ShannonDifferencesVectors$Distance <- abs(ShannonDifferencesVectors$Distance)
ShannonDifferencesVectors$Distance


BigData <- data.frame("Status" = BigDataJaccardsVectors1$Label, "PatientID" = BigdataYCVectors2$Label, 
                      "TimeDiff" = TimeDifferencesVectors$Distance, "SCORADDiff" = SCORADDifferencesVectors$Distance, 
                      "Jaccard Community similarity" = BigDataJaccardsVectors1$Distance, 
                      "YC Community similarity" = BigdataYCVectors1$Distance, 
                      "pHDifferences" = pHDifferencesVectors$Distance,
                      "TEWLDifferences" = TEWLDifferencesVectors$Distance,
                      "RichnessDifferences" = RichnessDifferencesVectors$Distance,
                      "ShannonDifference" = ShannonDifferencesVectors$Distance)

BigData <- BigData[grep("^Within", BigData$PatientID),]
BigData <- BigData[grep("^Within", BigData$Status),]

colourCount = length(unique(BigData$Status))
getPalette = colorRampPalette(brewer.pal(3, "Set1"))

BigData$Status <- droplevels(BigData$Status)
BigData$Status1 <- factor(BigData$Status, levels = c(levels = "Within AD.LS", "Within AD.NLS", "Within HC"),
  labels = c("AD.LS", "AD.NLS", "HC"))


InterindvidualVariation2a <- 
  ggplot(BigData, aes(x = BigData$PatientID, y = BigData$RichnessDifferences, 
  fill = BigData$Status1)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #facet_grid( ~ PatientID, scales = "free_x", space='free_x') + 
  scale_fill_manual(values = getPalette(colourCount)) + ggtitle("a") + labs(fill = NULL,  y = "Rich. Diff.") +
  theme(axis.title.x=element_blank())

InterindvidualVariation2b <- 
  ggplot(BigData, aes(x = BigData$PatientID, y = BigData$Jaccard.Community.similarity, 
  fill = BigData$Status1)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #facet_grid( ~ PatientID, scales = "free_x", space='free_x') + 
  scale_fill_manual(values = getPalette(colourCount)) + ggtitle("b") + labs(fill = NULL,  y = "Community mem.") +
  theme(axis.title.x=element_blank())

InterindvidualVariation2c <- 
  ggplot(BigData, aes(x = BigData$PatientID, y = BigData$YC.Community.similarity, 
  fill = BigData$Status1)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #facet_grid( ~ PatientID, scales = "free_x", space='free_x') + 
  scale_fill_manual(values = getPalette(colourCount)) + ggtitle("c") + labs(fill = NULL,  y = "Community sim.") +
  theme(axis.title.x=element_blank())

library(patchwork)
InverindividualVariationFigure <- InterindvidualVariation2a + InterindvidualVariation2b + InterindvidualVariation2c & theme(legend.position = "bottom")
InverindividualVariationFigure <- InverindividualVariationFigure + plot_layout(guides = "collect") + plot_layout(ncol = 3)
InverindividualVariationFigure


model <- lmer(BigData$RichnessDifferences ~ BigData$PatientID +
  (1|BigData$Status), data = Tape1MatrixData)
drop1(model, test = "Chisq")

model <- lmer(BigData$Jaccard.Community.similarity ~ BigData$PatientID +
  (1|BigData$Status), data = Tape1MatrixData)
drop1(model, test = "Chisq")

model <- lmer(BigData$YC.Community.similarity ~ BigData$PatientID +
  (1|BigData$Status), data = Tape1MatrixData)
drop1(model, test = "Chisq")











metadata$LS.NLS
summary(metadata$TapeNo)








##########Staph investigations#############

perfect_staph_hits <- c("SeqVar11000", "SeqVar11242", "SeqVar1315", "SeqVar13448",
"SeqVar13865", "SeqVar1408", "SeqVar146", "SeqVar147", "SeqVar15689",
"SeqVar1617", "SeqVar1619", "SeqVar162", "SeqVar18708", "SeqVar20728",
"SeqVar23", "SeqVar236", "SeqVar23779", "SeqVar2539", "SeqVar2653",
"SeqVar2659", "SeqVar285", "SeqVar29436", "SeqVar347", "SeqVar359",
"SeqVar368", "SeqVar3800", "SeqVar4013", "SeqVar5993", "SeqVar6042",
"SeqVar606", "SeqVar624", "SeqVar6282", "SeqVar719", "SeqVar7503",
"SeqVar771", "SeqVar7860", "SeqVar907", "SeqVar94")

row.names(Taxonomy)
a <- Taxonomy %>% filter(
  row.names(Taxonomy) %in% perfect_staph_hits
)

Taxaonomy$SeqID %in% perfect_staph_hits
a <- which(row.names(Taxonomy) %in% perfect_staph_hits)
rowSums(rel.data[, c(a)])



