library(recount3)
#Descargando todos los proyectos guardados en recount3 en una variable
human_projects <- available_projects()

#Obteniendo mi estudio de interés
rse_jana <- create_rse( subset(human_projects, project== "SRP113338" & project_type == "data_sources"))

# Modificando el contenido del assay a cuentas por lectura y no por nucléotido
assay(rse_jana, "raw_counts") <- compute_read_counts(rse_jana)

# Obteniendo más información de las muestras
rse_jana <- expand_sra_attributes(rse_jana) #Este comando no resulta por lo que se modifica la estructura de los atributos de la muestra
rse_jana$sra.sample_attributes[1:3]

# Se genera un archivo con los atributos de las muestras para evaluar el problema
namee <- "C:/Users/Acer/Documents/ibt_ccg/Bioinfo2/RNA-seq/colData_info_rseSRP113338.txt"
write.table(rse_jana$sra.sample_attributes, file=namee, row.names = TRUE)

# Dado que el problema está en las últimas líneas en donde falta una categoría se sustituye
# todo el contenido de la línea por una versión en done tenga " | NA "
rse_jana$sra.sample_attributes[7:8] <- "source_name;;normal liver tissue|tissue;;normal liver| NA"
rse_jana$sra.sample_attributes[9] <- "cell type;;primary hepatocytes|source_name;;primary human hepatocyte| NA"

# Volver a correr expand_sra_attributes(rse_jana)
rse_jana <- expand_sra_attributes(rse_jana)

# Dando el formato de R adecuado a cada uno de los vectores
rse_jana$sra_attribute.cell_type<- factor(rse_jana$sra_attribute.cell_type)
rse_jana$sra_attribute.day <- gsub("day", "", rse_jana$sra_attribute.day)
rse_jana$sra_attribute.day <- gsub("[^0-9.]+", "0", rse_jana$sra_attribute.day)
rse_jana$sra_attribute.day <- as.factor (rse_jana$sra_attribute.day)
rse_jana$sra_attribute.source_name<- factor(rse_jana$sra_attribute.source_name)

colData(rse_jana)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_jana)))
]

summary(as.data.frame(colData(rse_jana)[
  ,
  grepl("^sra_attribute.[cell_type|day|source_name]", colnames(colData(rse_jana)))
]))

# Evaluando el estado de la secuenciación del rse
rse_jana$assigned_gene_prop <- rse_jana$recount_qc.gene_fc_count_all.assigned / rse_jana$recount_qc.gene_fc_count_all.total
summary(rse_jana$assigned_gene_prop)

# Evaluando las muestras con la variable assigned_gene_prop
attributes <- colnames(colData(rse_jana))[grepl("^sra_attribute.[cell_type|day|source_name]",colnames(colData(rse_jana)))]
library(ggplot2)
hist(rse_jana$assigned_gene_prop)
data <- rse_jana$assigned_gene_prop

# Analizando cuantas muestras tienen una proporciones de genes menor a 0.4
table(rse_jana$assigned_gene_prop < 0.4)

#Evaluando la expresión de los genes
gene_means <- rowMeans(assay(rse_jana,"raw_counts"))
summary(gene_means)
rse_jana <- rse_jana [gene_means> 0.1, ] # Quedan 19708 genes después del filtrado

#Visualizando las medias de expresión
temp <- gene_means
temp <- as.matrix(temp)
View(temp)

# Normalizando los datos
library("edgeR")
# Crear un objeto de tipo DGElista, extrayendo todos los metadatos de los genes
dge <- DGEList(
  counts = assay(rse_jana, "raw_counts"),
  genes = rowData(rse_jana)
)
# Calcula los factores de normalización para cada una de las muestras
dge <- calcNormFactors(dge)

# Modelo estadístico

mod <- model.matrix(~ sra_attribute.cell_type,
                    data= colData(rse_jana)
)
#Visualización del modelo
library(ExploreModelMatrix)
vd <- VisualizeDesign(sampleData = colData(rse_jana),
                      designFormula = ~ sra_attribute.cell_type,
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
colnames(mod)

# Utilizando el modelo estadístico creado anteriormente y los datos de los genes
# de dge se crea un voom plot

library("limma")
vGene <- voom(dge, mod, plot = TRUE)

# Calculos bayesianos

eb_results <- eBayes(lmFit(vGene))

# Tabla con los resultados top de eb_results
de_results <- topTable(
  eb_results,
  number = nrow(rse_jana),
  sort.by = "none"
)

dim(de_results)
head(de_results)

# Análisis de parámetros estadísticos

table(de_results$adj.P.Val > 0.05)
write.table(de_results$adj.P.Val, file=ruta, row.names = TRUE)

#Más visualizaciones

# Generar los plotMA para cada uno de los modelos celulares evaluados
for (i in 2:5) {
  # Generar el plotMA y guardarlo como una imagen en la carpeta "plots"
  png(paste0("figuras/plot_", i, ".png"))
  plotMA(eb_results, coef=i)
  dev.off()
}

# Generar los volcanoplots para cada uno de los modelos celulares evaluados
for (i in 2:5) {
  # Generar el plotMA y guardarlo como una imagen en la carpeta "plots"
  png(paste0("figuras/volcanoplot_", i, ".png"))
  volcanoplot(eb_results, coef= i, highlight = 3, names = de_results$gene_name)
  dev.off()
}

#Heatmap en iSSEE

## Extraer valores de los genes de interés

exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]
df <- as.matrix(exprs_heatmap)

heatmap(exprs_heatmap)
iSEE(rse_jana)

# Extras
library("ggplot2")
ggplot(as.data.frame(colData(rse_jana)), aes(y = assigned_gene_prop, x = sra_attribute.cell_type)) +
  geom_boxplot() +
  theme_bw(base_size = 10) +
  ylab("Assigned Gene Prop") +
  xlab("Cell type")
