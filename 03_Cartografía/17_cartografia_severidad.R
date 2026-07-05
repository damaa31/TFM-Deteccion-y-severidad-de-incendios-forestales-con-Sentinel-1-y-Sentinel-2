# =============================================================================
# TFM SENTINEL-1 · CARTOGRAFIA DE SEVERIDAD (Objetivo 2)
# =============================================================================
# Mapas de severidad CBI predicha dentro del perimetro EFFIS (la severidad solo
# se define sobre area quemada; Key & Benson 2006, Miller & Thode 2007). Se usa
# EFFIS como mascara (delimitacion oficial e independiente) en vez de la
# deteccion SAR, para no propagar la comision del modelo de deteccion y evitar
# circularidad. El CBI bajo dentro del perimetro corresponde a zonas poco o nada
# afectadas (islas internas), coherente con la naturaleza continua del CBI.
#
# Para cada modelo (SAR, Optico): CBI continuo (0-3) y categorico (Baja/Mod/Alta).
# Modelo entrenado con los OTROS incendios (coherente con LOFOCV).
#
# COMPROBACIONES para la memoria:
#   1. Validacion contra parcelas CBI reales de Llamas (R2, RMSE) -> mapa vs campo.
#   2. Matriz de confusion categorica mapa vs clases reales.
#   3. Coincidencia espacial SAR vs Optico (kappa entre mapas).
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(sf); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  rds_opt_eco = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/optico_por_ecosistema_v3.rds",
  dir_s1      = "G:/TFM/_DATOS/EXPORT_S1_TFM",
  dir_s2      = "G:/TFM/_DATOS/EXPORT_S2_TFM",
  shp_effis   = "G:/TFM/_CARTOGRAFIA/PERIMETROS_INCENDIOS/Perimetros_CyL.shp",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_CARTOGRAFIA_SEVERIDAD/Porto"
)
INCENDIO <- "Porto"
N_ARBOLES <- 500
SEED <- 42
ETIQUETAS <- c("Baja","Moderada","Alta")
dir.create(RUTAS$dir_salida, showWarnings=FALSE, recursive=TRUE)

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
opt_eco <- readRDS(RUTAS$rds_opt_eco)
dataset <- rds$dataset
bandas_sar <- coln$bandas_depuradas
bandas_opt <- opt_eco$bandas_opt_dep

cat("=============================================================\n")
cat(sprintf("CARTOGRAFIA DE SEVERIDAD · %s (mascara EFFIS)\n", INCENDIO))
cat("=============================================================\n")

effis <- st_read(RUTAS$shp_effis, quiet=TRUE)
if (st_crs(effis)$epsg != 25830) effis <- st_transform(effis, 25830)
campo <- NULL
for (c in c("COMMUNE","Incendio","NOMBRE","name")) if (c %in% names(effis)) {campo<-c; break}
per_inc <- st_union(effis[grepl(INCENDIO, effis[[campo]], ignore.case=TRUE), ])

clasificar_cbi <- function(x) {
  cut(x, breaks=c(-0.01,1.0,2.0,3.01), labels=ETIQUETAS)
}
kappa_cohen <- function(a, b) {
  a<-factor(a,levels=ETIQUETAS); b<-factor(b,levels=ETIQUETAS)
  cm<-table(a,b); n<-sum(cm); acc<-sum(diag(cm))/n
  pe<-sum(rowSums(cm)*colSums(cm))/n^2
  (acc-pe)/(1-pe)
}
buscar_tif <- function(d,inc){t<-list.files(d,pattern="\\.tif$",full.names=TRUE)
m<-t[grepl(paste0(inc,"\\.tif"),basename(t))]; if(length(m)==0) NA_character_ else m[1]}

# Modelo + mapa, devuelve el raster continuo
cartografiar <- function(dir_tif, bandas, etiqueta) {
  cat(sprintf("\n--- Modelo %s ---\n", etiqueta))
  tif <- buscar_tif(dir_tif, INCENDIO)
  if (is.na(tif)) { cat("  Sin TIF.\n"); return(NULL) }
  train <- dataset[dataset$Incendio != INCENDIO, ]
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  set.seed(SEED)
  rf <- ranger(form, data=train, num.trees=N_ARBOLES, seed=SEED)
  r <- rast(tif); r <- r[[intersect(bandas, names(r))]]
  pf <- function(model,data){p<-predict(model,data,num.threads=1)$predictions; pmin(pmax(p,0),3)}
  cat("  Prediciendo...\n")
  cbi <- predict(r, rf, fun=pf, na.rm=TRUE); names(cbi)<-"CBI"
  cbi <- mask(crop(cbi, vect(per_inc)), vect(per_inc))
  cat_map <- classify(cbi, matrix(c(-0.01,1,1, 1,2,2, 2,3.01,3),ncol=3,byrow=TRUE))
  writeRaster(cbi, file.path(RUTAS$dir_salida, sprintf("sev_CBI_%s_%s.tif",etiqueta,INCENDIO)), overwrite=TRUE)
  writeRaster(cat_map, file.path(RUTAS$dir_salida, sprintf("sev_clase_%s_%s.tif",etiqueta,INCENDIO)), overwrite=TRUE)
  v<-values(cbi); v<-v[!is.na(v)]
  cat(sprintf("  CBI predicho: media=%.2f rango=[%.2f,%.2f]\n", mean(v),min(v),max(v)))
  list(modelo=rf, cbi=cbi)
}

res_sar <- cartografiar(RUTAS$dir_s1, bandas_sar, "SAR")
res_opt <- cartografiar(RUTAS$dir_s2, bandas_opt, "Optico")

# =============================================================================
# COMPROBACIONES PARA LA MEMORIA
# =============================================================================
cat("\n=============================================================\n")
cat("COMPROBACIONES (validacion del mapa)\n")
cat("=============================================================\n")

# Parcelas reales de Llamas (solo campo, con CBI observado)
parc <- dataset[dataset$Incendio==INCENDIO & dataset$origen=="campo", ]
parc_v <- vect(parc, geom=c("x_utm","y_utm"), crs="EPSG:25830")
cat(sprintf("Parcelas de campo en %s: %d\n", INCENDIO, nrow(parc)))

validar_mapa <- function(res, etiqueta) {
  if (is.null(res)) return(NULL)
  pred <- terra::extract(res$cbi, parc_v)[,2]
  obs  <- parc$CBI_total
  ok <- !is.na(pred)
  pred<-pred[ok]; obs<-obs[ok]
  # Continuo
  r2 <- 1 - sum((obs-pred)^2)/sum((obs-mean(obs))^2)
  rmse <- sqrt(mean((obs-pred)^2))
  cat(sprintf("\n--- %s vs parcelas reales (n=%d) ---\n", etiqueta, length(obs)))
  cat(sprintf("  Continuo: R2=%.3f | RMSE=%.3f | r=%.3f\n", r2, rmse, cor(obs,pred)))
  # Categorico
  cm <- table(Obs=clasificar_cbi(obs), Pred=clasificar_cbi(pred))
  acc <- sum(diag(cm))/sum(cm)
  cat(sprintf("  Categorico: Accuracy=%.3f | Kappa=%.3f\n", acc,
              kappa_cohen(clasificar_cbi(obs), clasificar_cbi(pred))))
  cat("  Matriz de confusion (obs filas / pred col):\n"); print(cm)
  data.frame(modelo=etiqueta, n=length(obs), r2=round(r2,3),
             rmse=round(rmse,3), accuracy=round(acc,3))
}

v_sar <- validar_mapa(res_sar, "SAR")
v_opt <- validar_mapa(res_opt, "Optico")

# Coincidencia espacial SAR vs Optico (mapas categoricos)
if (!is.null(res_sar) && !is.null(res_opt)) {
  cat("\n--- Coincidencia espacial SAR vs Optico (mapa categorico) ---\n")
  cs <- values(classify(res_sar$cbi, matrix(c(-0.01,1,1,1,2,2,2,3.01,3),ncol=3,byrow=TRUE)))
  co <- values(classify(res_opt$cbi, matrix(c(-0.01,1,1,1,2,2,2,3.01,3),ncol=3,byrow=TRUE)))
  ok <- !is.na(cs) & !is.na(co)
  coincidencia <- mean(cs[ok]==co[ok])
  ks <- kappa_cohen(ETIQUETAS[cs[ok]], ETIQUETAS[co[ok]])
  cat(sprintf("  Coincidencia de clase: %.1f%% | Kappa entre mapas: %.3f\n",
              100*coincidencia, ks))
  # Diferencia continua
  dif <- res_opt$cbi - res_sar$cbi
  writeRaster(dif, file.path(RUTAS$dir_salida, sprintf("sev_dif_OptMenosSAR_%s.tif",INCENDIO)), overwrite=TRUE)
  vd<-values(dif); vd<-vd[!is.na(vd)]
  cat(sprintf("  Diferencia media Optico-SAR: %.2f CBI\n", mean(vd)))
}

cat("\n--- Resumen validacion ---\n")
print(rbind(v_sar, v_opt), row.names=FALSE)
cat(sprintf("\nSalidas en: %s\n", RUTAS$dir_salida))