# =============================================================================
# TFM SENTINEL-1 · CARTOGRAFIA DE SEVERIDAD — SOLO FUSION (Objetivo 2)
# =============================================================================
# Genera UNICAMENTE el TIF de fusion (SAR + Optico) de CBI continuo y categorico
# en la carpeta de cada incendio. No toca los TIF de SAR/Optico ya existentes.
# Modelo entrenado con los OTROS incendios (LOFOCV), mascara EFFIS.
# La fusion combina bandas de DOS TIF distintos (S1 y S2), por eso se cargan y
# se apilan ambos rasters antes de predecir.
# Reporta validacion contra parcelas de campo (R2, RMSE, kappa categorico).
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(sf); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_optfus  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/optico_fusion_severidad_v3.rds",
  dir_s1      = "G:/TFM/_DATOS/EXPORT_S1_TFM",
  dir_s2      = "G:/TFM/_DATOS/EXPORT_S2_TFM",
  shp_effis   = "G:/TFM/_CARTOGRAFIA/PERIMETROS_INCENDIOS/Perimetros_CyL.shp",
  dir_base    = "G:/TFM/_DATOS/_ANALISIS/_CARTOGRAFIA_SEVERIDAD"  # contiene 8 subcarpetas
)

# Los 8 incendios = nombres de las subcarpetas
INCENDIOS <- c("Resoba","Llamas","Molezuelas","Orallo","Fasgar","Canalejas","Barniedo","Porto")
N_ARBOLES <- 500
SEED <- 42
ETIQUETAS <- c("Baja","Moderada","Alta")

# --- Datos y bandas de fusion (8 SAR + 6 opticas = 14, ya guardadas en el rds) ---
rds    <- readRDS(RUTAS$rds_dataset)
optfus <- readRDS(RUTAS$rds_optfus)
dataset       <- rds$dataset
bandas_fusion <- optfus$bandas_fusion   # 14 bandas (SAR depuradas + opticas depuradas)

# --- Perimetros EFFIS ---
effis <- st_read(RUTAS$shp_effis, quiet=TRUE)
if (st_crs(effis)$epsg != 25830) effis <- st_transform(effis, 25830)
campo <- NULL
for (c in c("COMMUNE","Incendio","NOMBRE","name")) if (c %in% names(effis)) {campo<-c; break}

clasificar_cbi <- function(x) cut(x, breaks=c(-0.01,1.0,2.0,3.01), labels=ETIQUETAS)
kappa_cohen <- function(a, b) {
  a<-factor(a,levels=ETIQUETAS); b<-factor(b,levels=ETIQUETAS)
  cm<-table(a,b); n<-sum(cm); acc<-sum(diag(cm))/n
  pe<-sum(rowSums(cm)*colSums(cm))/n^2
  (acc-pe)/(1-pe)
}
buscar_tif <- function(d,inc){t<-list.files(d,pattern="\\.tif$",full.names=TRUE)
m<-t[grepl(paste0(inc,"\\.tif"),basename(t))]; if(length(m)==0) NA_character_ else m[1]}

resumen <- data.frame()

for (INCENDIO in INCENDIOS) {
  cat("\n=============================================================\n")
  cat(sprintf("FUSION · %s\n", INCENDIO))
  cat("=============================================================\n")
  
  dir_salida <- file.path(RUTAS$dir_base, INCENDIO)
  if (!dir.exists(dir_salida)) { cat("  Carpeta no encontrada, salto.\n"); next }
  
  # Perimetro del incendio
  per_inc <- st_union(effis[grepl(INCENDIO, effis[[campo]], ignore.case=TRUE), ])
  
  # --- Cargar y APILAR los dos TIF (S1 + S2) ---
  tif_s1 <- buscar_tif(RUTAS$dir_s1, INCENDIO)
  tif_s2 <- buscar_tif(RUTAS$dir_s2, INCENDIO)
  if (is.na(tif_s1) || is.na(tif_s2)) { cat("  Falta TIF S1 o S2, salto.\n"); next }
  
  r1 <- rast(tif_s1)
  r2 <- rast(tif_s2)
  # Alinear S2 a la malla de S1 si difieren (mismo CRS/resolucion/extent)
  if (!compareGeom(r1, r2, stopOnError=FALSE)) {
    r2 <- resample(r2, r1, method="bilinear")
  }
  r <- c(r1, r2)                                  # apila todas las bandas
  # Comprobacion: deben estar las 14 ANTES de seleccionar
  faltan <- setdiff(bandas_fusion, names(r))
  if (length(faltan) > 0) { cat("  FALTAN bandas:", paste(faltan,collapse=", "), "- salto.\n"); next }
  r <- r[[bandas_fusion]]                         # selecciona EN EL ORDEN EXACTO del modelo
  
  # --- Entrenar modelo fusion con los OTROS incendios (LOFOCV) ---
  train <- dataset[dataset$Incendio != INCENDIO, ]
  form  <- as.formula(paste("CBI_total ~", paste(bandas_fusion, collapse=" + ")))
  set.seed(SEED)
  rf <- ranger(form, data=train, num.trees=N_ARBOLES, seed=SEED)
  
  # --- Predecir mapa (escribe por bloques a archivo: pico de memoria bajo) ---
  pf <- function(model, data, ...){ p<-predict(model,data,num.threads=1)$predictions; pmin(pmax(p,0),3) }
  cat("  Prediciendo fusion...\n")
  tmp_pred <- file.path(dir_salida, sprintf("tmp_pred_%s.tif", INCENDIO))
  cbi <- predict(r, rf, fun=pf, na.rm=TRUE,
                 filename=tmp_pred, overwrite=TRUE, wopt=list(steps=40))
  names(cbi) <- "CBI"
  cbi <- mask(crop(cbi, vect(per_inc)), vect(per_inc))
  cat_map <- classify(cbi, matrix(c(-0.01,1,1, 1,2,2, 2,3.01,3), ncol=3, byrow=TRUE))
  
  # --- Guardar SOLO los TIF de fusion ---
  writeRaster(cbi,     file.path(dir_salida, sprintf("sev_CBI_Fusion_%s.tif",  INCENDIO)), overwrite=TRUE)
  writeRaster(cat_map, file.path(dir_salida, sprintf("sev_clase_Fusion_%s.tif",INCENDIO)), overwrite=TRUE)
  v<-values(cbi); v<-v[!is.na(v)]
  cat(sprintf("  CBI predicho (fusion): media=%.2f rango=[%.2f,%.2f]\n", mean(v),min(v),max(v)))
  
  # --- Validacion contra parcelas de campo ---
  parc <- dataset[dataset$Incendio==INCENDIO & dataset$origen=="campo", ]
  parc_v <- vect(parc, geom=c("x_utm","y_utm"), crs="EPSG:25830")
  pred <- terra::extract(cbi, parc_v)[,2]
  obs  <- parc$CBI_total
  ok <- !is.na(pred); pred<-pred[ok]; obs<-obs[ok]
  r2v  <- 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
  rmse <- sqrt(mean((obs-pred)^2))
  cm   <- table(Obs=clasificar_cbi(obs), Pred=clasificar_cbi(pred))
  acc  <- sum(diag(cm))/sum(cm)
  kap  <- kappa_cohen(clasificar_cbi(obs), clasificar_cbi(pred))
  cat(sprintf("  Fusion vs campo (n=%d): R2=%.3f | RMSE=%.3f | r=%.3f | Acc=%.3f | Kappa=%.3f\n",
              length(obs), r2v, rmse, cor(obs,pred), acc, kap))
  
  resumen <- rbind(resumen, data.frame(
    Incendio=INCENDIO, n=length(obs), R2=round(r2v,3),
    RMSE=round(rmse,3), r=round(cor(obs,pred),3),
    Acc_cat=round(acc,3), Kappa_cat=round(kap,3)
  ))
  
  # Liberar memoria antes del siguiente incendio
  rm(r, r1, r2, cbi, cat_map, rf); gc()
  if (file.exists(tmp_pred)) file.remove(tmp_pred)
}

cat("\n=============================================================\n")
cat("RESUMEN FUSION (todos los incendios)\n")
cat("=============================================================\n")
print(resumen, row.names=FALSE)
cat(sprintf("\nTIF de fusion guardados en cada subcarpeta de: %s\n", RUTAS$dir_base))