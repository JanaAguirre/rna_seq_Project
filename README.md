---
title: "Reporte_proyecto"
author: "Jana A. Castañeda"
date: "2024-02-13"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Proyecto RNA-seq 2024
### Jana Carolina Aguirre Castañeda.
 12 - febrero - 2024

Sobre el proyecto
Este proyecto busca poner en práctica los conocimientos adquiridos durante el módulo Análisis de datos de secuenciación masiva del curso Bioinformática y estadística II.

## Estudio seleccionado: SRP113338
Los datos utilizados para este proyecto fueron obtenidos del estudio ** Whole transcriptome analysis of iPSC-derived hepatic organoids cultures **.
*Abstract*: "Overall design: Whole transcriptome analysis was performed on iPSCs, primary human hepatocytes (PHH), HO1 (primary hepatic organoids) and human liver tissue. For this analysis, over 100 organoids that were generated from a single donor were analyzed. Next, whole transcriptome analysis was performed to characterize HO2s (secondary organoids) and their relationship to iPSCs, PHHs and human liver tissue. This analysis utilized 1,000 HO2s that were prepared from a single donor. A heatmap, which focused upon 2458 genes that were involved in stem cell differentiation and liver development, was prepared to visualize the relative relationship between the transcriptomes of iPSCs, PHHs, human liver tissue and HO1. The dendrogram indicates that HO1 is more closely related to liver tissue than to PHHs; this is because cholangiocytes are present in HO1 and liver, but are absent (or rare) in PHHs. Examination of the genes within this heatmap indicated that there was a lack of expression of pluripotent genes in HO1; and that HO1 had an increased level of expression of cellular markers for liver lineages and liver development, which included mature hepatocytes, bile duct cells and NOTCH signaling components. A cluster analysis of 1510 genes indicates that there were major changes in the transcriptome after the cells that were dissociated from an HO1 were transferred from the growth medium into the differentiation medium. As with the HO1s, the gene expression profile of HO2s resembled that of PHHs and human liver. Of importance, the level of expression of mRNAs that are of importance for various human liver functions (ADH4, CYP3A4, TTR, GSRA1, TDO2, GSTA1 and FAH) were markedly increased in HO2s."

# 1 Librerías a utilizar.

```{r}
library("recount3")
library(ExploreModelMatrix)
library("edgeR")
library("limma")
library("ggplot2")
```
# 2 Datos a utilizar. 
Descargar los datos de RNA-seq desde el paquete recount3.
```{r}
#Descargando todos los proyectos guardados en recount3 en una variable
human_projects <- available_projects()

#Obteniendo el estudio de interés
rse_jana <- create_rse( subset(human_projects, project== "SRP113338" & project_type == "data_sources"))
# Modificando el contenido del assay a cuentas por lectura y no por nucléotido
assay(rse_jana, "raw_counts") <- compute_read_counts(rse_jana)
```
# 3 Formateo de los datos. 
Podemos observar los atributos que contiene nuestro objeto ahora.
```{r}
rse_jana$sra.sample_attributes[1:3]
```
Hasta esta punto parece que las dimensiones de los vectores son correctas, por
lo que se podría proseguir a expandir los atributos de las muestras. 
```{r}
rse_jana <- expand_sra_attributes(rse_jana)
```
Al correr ``expand_sra_attributes(rse_jana)`` se genera una advertencia y la
acción no se ejecuta, por lo tanto se prosigue a averiguar el problema generando
un archivo con contenga todos los elementos del vector "rse_jana$sra.sample_attributes".
```{r}
# Se genera un archivo con los atributos de las muestras para evaluar el problema
namee <- "C:/Users/Acer/Documents/ibt_ccg/Bioinfo2/RNA-seq/colData_info_rseSRP113338.txt"
write.table(rse_jana$sra.sample_attributes, file=namee, row.names = TRUE)

```
Al abrir el archivo
```{r}
# Lee el archivo de texto
contenido <- readLines("colData_info_rseSRP113338.txt")
print(contenido)
```
Se puede apreciar que las muestras 7,8 y 9 solo tienen dos columnas. Con esta información se pasa a sustituir ese contenido por líneas que tengan una tercera
columna igual a " | NA ".
```{r}
# Dado que el problema está en las últimas líneas en donde falta una categoría se sustituye
# todo el contenido de la línea por una versión en done tenga " | NA "
rse_jana$sra.sample_attributes[7:8] <- "source_name;;normal liver tissue|tissue;;normal liver| NA"
rse_jana$sra.sample_attributes[9] <- "cell type;;primary hepatocytes|source_name;;primary human hepatocyte| NA"
```
Ahora la función ``expand_sra_attributes(rse_jana)`` funciona bien
```{r}
# Volver a correr expand_sra_attributes(rse_jana)
rse_jana <- expand_sra_attributes(rse_jana)
```
## 3.1 Asignación de tipo de vector apropiado.

En este estudio las únicas variables estudiadas son los distintos tipos de organoides derivados de iPSC de un mismo individuo. Por lo tanto, el único factor a considerar son las diferencias transcripcionales entre cada uno de los grupos. Para esto, en la columna sra_attribute_cell_type debe ser considerada como un factor y no una columna de caracteres. Con el fin de que análisis posteriores cada uno de los elementos corresponda a un "nivel" o "level" dentro del factor. 

```{r}
# Dando el formato de R adecuado a cada uno de los vectores
rse_jana$sra_attribute.cell_type<- factor(rse_jana$sra_attribute.cell_type)
rse_jana$sra_attribute.day <- gsub("day", "", rse_jana$sra_attribute.day)
rse_jana$sra_attribute.day <- gsub("[^0-9.]+", "0", rse_jana$sra_attribute.day)
rse_jana$sra_attribute.day <- as.factor (rse_jana$sra_attribute.day)
rse_jana$sra_attribute.source_name<- factor(rse_jana$sra_attribute.source_name)
```
También se modifica el elemento "sra_attribute.day" por motivos de practicidad en caso de que se quisiera ocupar el elemento. 

# 4. Filtrado de datos
Calcularemos la proporción de lecturas asignadas a genes
# Evaluando el estado de la secuenciación del rse

```{r}
rse_jana$assigned_gene_prop <- rse_jana$recount_qc.gene_fc_count_all.assigned / rse_jana$recount_qc.gene_fc_count_all.total
summary(rse_jana$assigned_gene_prop)
```
```{r rse_jana$assigned_gene_prop, echo= FALSE}
# Evaluando las muestras con la variable assigned_gene_prop

hist(rse_jana$assigned_gene_prop)
```
```{r}
table(rse_jana$assigned_gene_prop < 0.4)
```
La cantidad de muestras es muy pequeña (n=9), de ahí que las frecuencias del histograma no sobrepasan 3. Por otro lado, podemos observar que practicamente a la mitad de las muestras no se les asignó un alto porcentage de los genes totales detectados. Y dado que todos los organoides deberían tener un background genético y de expresión muy similar, podemos asumir que estas 4 muestras podrían no haber tenido la mejor calidad. Sin embargo, dado la n de la muestra y a que no hay replicas para practicamente ningún tipo de organoide se mantendrán las muestras con discresión.  

### 4.1 Filtrado de los genes
Dado que hay una gran cantidad de genes cuya expresión es 0 y alterarían los análisis estadísticos se filtran.
```{r}
#Evaluando la expresión de los genes
gene_means <- rowMeans(assay(rse_jana,"raw_counts"))
rse_jana <- rse_jana [gene_means> 0.1, ] 
summary(gene_means) # Quedan 19708 genes después del filtrado

#Visualizando las medias de expresión
temp <- gene_means
temp <- as.matrix(temp)
View(temp)
```
# 5. Normalización de datos

Para reducir la incidencia de falsos positivos normalizaremos los datos asumiendo que la mayoría de los genes no se están expresando diferencialmente.
```{r}
library("edgeR")
# Crear un objeto de tipo DGElista, extrayendo todos los metadatos de los genes
dge <- DGEList(
  counts = assay(rse_jana, "raw_counts"),
  genes = rowData(rse_jana)
)
# Calcula los factores de normalización para cada una de las muestras
dge <- calcNormFactors(dge)
``` 
# 6. Expresión diferencial
# 6.1 Modelo estadístico
Para este estudio los atributos que se generaron consistían solo de columnas categóricas en donde se mencionaba el tipo de organoide y en algunos casos se registraron días durante los que se dejaron crecer. Por lo tanto, en este modelo estadístico se plantea que "Y" solo depende de sra_attribute.cell_type. 
```{r}
# Modelo estadístico

mod <- model.matrix(~ sra_attribute.cell_type,
                    data= colData(rse_jana)
)
```

El modelo queda más claro al visualizarlo.
```{r echo=FALSE}
#Visualización del modelo
library(ExploreModelMatrix)
vd <- VisualizeDesign(sampleData = colData(rse_jana),
                      designFormula = ~ sra_attribute.cell_type,
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
```
Aquí podemos observar que la variable de referencia son las muestras de hepatocitos, normales. Por lo que el modelo queda como el coeficiente de las muestras de hepatocitos - el coeficiente de alguno de los grupos de organoides evaluados. 

```{r}
#Visualizando las columnas del modelo 
colnames(mod)
```
Utilizando el modelo estadístico creado anteriormente y los features de los genes contenidos en la variable dge se crea un voom plot.
```{r}

library("limma")
vGene <- voom(dge, mod, plot = TRUE)

```
La desviación estándar de las muestras que tuvieron log2 mayores a 5 es grande, por esta razón en lugar de observar una decaída de la función como de costumbra se ve una línea ascendente. Sin embargo, la primera sección de los genes tienen el comportamiento de "J" esperado. 

# 7. Otras visualizaciones 

## plotMAs

Se generan plotMAs para cada uno de los cell.types para ver cómo se comporta el log-fold-change en cada grupo en función del average-log-expression.

```{r eval=FALSE, include=FALSE}
# Generar los plotMA para cada uno de los modelos celulares evaluados
for (i in 2:5) {
 Generar el plotMA y guardarlo como una imagen en la carpeta "plots"
 png(paste0("figuras/plot_", i, ".png"))
plotMA(eb_results, coef=i)
dev.off()}
```
![plots_RNAseq_all](https://github.com/JanaAguirre/rna_seq_Project/assets/158189687/f234dcc7-1ef3-40ac-bf25-1a1ae1a83923)

## volcanoplots

Con la misma lógica se pueden generar los volcanoplots para cada uno de los grupos. 
```{r eval=FALSE, include=FALSE}
# Generar los volcanoplots para cada uno de los modelos celulares evaluados
for (i in 2:5) {
   Generar el plotMA y guardarlo como una imagen en la carpeta "plots"
  png(paste0("figuras/volcanoplot_", i, ".png"))
  volcanoplot(eb_results, coef= i, highlight = 3, names = de_results$gene_name)
  dev.off()}
``` 
![volcanoplots_RNAseq_all](https://github.com/JanaAguirre/rna_seq_Project/assets/158189687/ee3d69bf-fa38-4bc3-bc97-50ac0c752a6f)

# 8. Evaluación estadística 
Para obtener los valores de pvalue así como otros parámetros estadísticos de intéres, se utilizó la función ```eBayes```.

```{r}
# Cálculos bayesianos

eb_results <- eBayes(lmFit(vGene))

# Tabla con los resultados top de eb_results
de_results <- topTable(
  eb_results,
  number = nrow(rse_jana),
  sort.by = "none"
)

dim(de_results)
```
Para que cada uno de nuestros "cell types" aparezcan en el "de_results" no específicamos el coeficiente. Con estos nuevos datos sobre la expresión diferencial del modelo, se puede evaluar cuantos de los genes tuvieron un pvalue menor a 0.05.
```{r}
# Análisis de parámetros estadísticos

table(de_results$adj.P.Val < 0.05)
# Opcionalmente se puede generar un archivo con esta información
ruta <- "C:/Users/Acer/Documents/ibt_ccg/Bioinfo2/RNA-seq/de_RESULTS_padjustvalue.txt"
write.table(de_results$adj.P.Val, file=ruta, row.names = TRUE)

```
Esta parte resulta interesante, dado que es convención pensar en la necesidad de pvalues menores a 0.05 para tener resultados estadísticos correctos y confiables. Sin embargo, en este estudio en particular se debe hacer notar que todos los organoides derivaban de un mismo paciente y el objetivo del desarrollo de este tipo de modelos de cultivo tridimensionales es reflejar el contexto transcriptómico, metabólico y de señalización del paciente de origen. Un pvalue menor a 0.05 indicaría que hay una probabilidad muy baja de que las medias entre los dos grupos sean iguales, pero en este contexto sí se esperan que sean iguales. Por lo tanto, los valores de pvalue elevados nos hablan de un buen reflejo transcriptómico de los modelos.  

## Conclusiones

El dataset con ID SRP113338 no solo da paso a la evaluación de los genes con expresión diferencial entre cada uno de los organoides generados (iPSCs, primary hepatic organoids, secondary hepatic organoids, primary hepatocytes), si no también ha permitido evidenciar lo robusto del método de modelado. Para evaluar los genes que tienen expresión particular en cada uno de estos grupos en estudios futuros, se propone utilizar a los 50 genes con log-fold-change más grandes o más bajos en vez de utilizar los valores de pvalue.
