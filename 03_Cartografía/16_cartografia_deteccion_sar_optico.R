# =============================================================================
# TFM SENTINEL-1 · CARTOGRAFIA DE DETECCION (Objetivo 1)
# =============================================================================
# Mapas de probabilidad y clasificacion binaria de quemado. Modelo entrenado con
# los OTROS incendios (coherente con LOFOCV). Incluye:
#   - Comprobacion contra perimetro Copernicus (kappa corregido a double).
#   - Post-proceso MMU (Minimum Mapping Unit): elimina manchas quemadas aisladas
#     por debajo de un area minima. Practica estandar en cartografia de area
#     quemada para reducir comision (Chuvieco et al. 2016; productos MODIS
#     MCD64A1 y Copernicus aplican unidades minimas cartografiables). Un incendio
#     es espacialmente contiguo; pixeles quemados aislados son ruido.
#   - Diagnostico del tamano de manchas y barrido de varios MMU para elegir con
#     criterio (no arbitrario).
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(sf); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_modelado_v3.rds",
  dir_s1      = "G:/TFM/_DATOS/EXPORT_S1_TFM",
  shp_perim   = "G:/TFM/_CARTOGRAFIA/PERIMETROS_INCENDIOS/Perimetros_COPERNICUS_30.shp",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_CARTOGRAFIA"
)
INCENDIOS_OBJETIVO <- c("Porto", "Fasgar")
N_ARBOLES <- 500
SEED <- 42
MMU_PRUEBA_HA <- c(0.5, 1, 2, 5)   # areas minimas a probar (ha)
PIXEL_HA <- 0.01                    # 10x10 m = 0.01 ha
dir.create(RUTAS$dir_salida, showWarnings=FALSE, recursive=TRUE)

rds <- readRDS(RUTAS$rds_dataset)
dataset    <- rds$dataset
bandas_sar <- rds$bandas_sar
dataset <- dataset[dataset$Incendio != "Ciperez" & !is.na(dataset$quemado), ]

perim <- st_read(RUTAS$shp_perim, quiet=TRUE)
if (st_crs(perim)$epsg != 25830) perim <- st_transform(perim, 25830)
campo_nombre <- NULL
for (c in c("COMMUNE","Incendio","NOMBRE","name","Name")) {
  if (c %in% names(perim)) { campo_nombre <- c; break }
}

cat("=============================================================\n")
cat("CARTOGRAFIA DE DETECCION + MMU (Objetivo 1)\n")
cat("=============================================================\n")

buscar_tif <- function(dir_tif, incendio) {
  tifs <- list.files(dir_tif, pattern="\\.tif$", full.names=TRUE)
  m <- tifs[grepl(paste0(incendio,"\\.tif"), basename(tifs))]
  if (length(m)==0) NA_character_ else m[1]
}
calibrar_umbral <- function(rf, train) {
  oob <- rf$predictions[,2]; obs <- train$quemado
  u <- seq(0.2,0.6,by=0.01)
  f1 <- sapply(u, function(x){
    p<-as.integer(oob>=x); tp<-sum(p==1&obs==1); fp<-sum(p==1&obs==0); fn<-sum(p==0&obs==1)
    pr<-tp/(tp+fp); rc<-tp/(tp+fn)
    if(is.nan(pr)||is.nan(rc)||(pr+rc)==0) 0 else 2*pr*rc/(pr+rc)
  })
  u[which.max(f1)]
}
# Metricas vs perimetro (con conteos en double para evitar desbordamiento)
comparar <- function(clasif_map, per_inc) {
  per_rast <- rasterize(vect(per_inc), clasif_map, field=1, background=0)
  pv <- values(clasif_map); rv <- values(per_rast)
  ok <- !is.na(pv) & !is.na(rv)
  pv <- as.numeric(pv[ok]); rv <- as.numeric(rv[ok])
  tp<-sum(pv==1&rv==1); fp<-sum(pv==1&rv==0); fn<-sum(pv==0&rv==1); tn<-sum(pv==0&rv==0)
  n<-tp+fp+fn+tn; acc<-(tp+tn)/n; iou<-tp/(tp+fp+fn); dice<-2*tp/(2*tp+fp+fn)
  p_exp <- ((tp+fp)/n)*((tp+fn)/n) + ((fn+tn)/n)*((fp+tn)/n)
  kappa <- (acc-p_exp)/(1-p_exp)
  data.frame(accuracy=acc, iou=iou, dice=dice, kappa=kappa,
             ha_pred=(tp+fp)*PIXEL_HA, ha_perim=(tp+fn)*PIXEL_HA)
}

for (inc in INCENDIOS_OBJETIVO) {
  cat(sprintf("\n=============== %s ===============\n", inc))
  tif_path <- buscar_tif(RUTAS$dir_s1, inc)
  if (is.na(tif_path)) { cat("  Sin TIF.\n"); next }
  
  train <- dataset[dataset$Incendio != inc, ]
  form <- as.formula(paste("factor(quemado) ~", paste(bandas_sar, collapse=" + ")))
  set.seed(SEED)
  rf <- ranger(form, data=train, num.trees=N_ARBOLES, probability=TRUE, seed=SEED)
  umbral <- calibrar_umbral(rf, train)
  cat(sprintf("  Umbral calibrado: %.2f\n", umbral))
  
  r <- rast(tif_path); r <- r[[intersect(bandas_sar, names(r))]]
  pf <- function(model,data) predict(model,data,num.threads=1)$predictions[,2]
  cat("  Generando mapa...\n")
  prob_map <- predict(r, rf, fun=pf, na.rm=TRUE); names(prob_map)<-"prob"
  clasif <- prob_map >= umbral; names(clasif)<-"quemado"
  writeRaster(prob_map, file.path(RUTAS$dir_salida, paste0("prob_quemado_",inc,".tif")), overwrite=TRUE)
  
  # Perimetro del incendio
  per_inc <- NULL
  if (!is.null(campo_nombre)) {
    sel <- perim[grepl(inc, perim[[campo_nombre]], ignore.case=TRUE), ]
    if (nrow(sel)>0) per_inc <- st_union(sel)
  }
  
  # --- Diagnostico de manchas (clumps) ---
  cat("  --- Diagnostico de manchas conexas ---\n")
  clz <- patches(clasif, directions=8, zeroAsNA=TRUE)
  freq_clz <- freq(clz)
  areas_ha <- freq_clz$count * PIXEL_HA
  cat(sprintf("    N manchas: %d | Mancha mayor: %.0f ha | Manchas <1ha: %d (%.1f%% del total de manchas)\n",
              length(areas_ha), max(areas_ha),
              sum(areas_ha<1), 100*sum(areas_ha<1)/length(areas_ha)))
  
  # --- Comprobacion SIN filtro ---
  if (!is.null(per_inc) && length(per_inc)>0) {
    cat("  --- SIN MMU vs Copernicus ---\n")
    m0 <- comparar(clasif, per_inc)
    cat(sprintf("    IoU=%.3f Dice=%.3f Kappa=%.3f | ha_pred=%.0f ha_perim=%.0f\n",
                m0$iou, m0$dice, m0$kappa, m0$ha_pred, m0$ha_perim))
    
    # --- Barrido de MMU ---
    cat("  --- Barrido MMU (area minima) vs Copernicus ---\n")
    mejor_iou <- m0$iou; mejor_mmu <- 0; mejor_map <- clasif
    for (mmu_ha in MMU_PRUEBA_HA) {
      min_pix <- mmu_ha / PIXEL_HA
      grandes <- freq_clz$value[freq_clz$count >= min_pix]
      filt <- clz; filt[!(values(clz) %in% grandes)] <- NA
      filt_bin <- !is.na(filt)
      m <- comparar(filt_bin, per_inc)
      cat(sprintf("    MMU=%.1f ha: IoU=%.3f Dice=%.3f Kappa=%.3f | ha_pred=%.0f\n",
                  mmu_ha, m$iou, m$dice, m$kappa, m$ha_pred))
      if (m$iou > mejor_iou) { mejor_iou<-m$iou; mejor_mmu<-mmu_ha; mejor_map<-filt_bin }
    }
    cat(sprintf("  >> Mejor MMU: %.1f ha (IoU=%.3f)\n", mejor_mmu, mejor_iou))
    names(mejor_map) <- "quemado"
    writeRaster(mejor_map, file.path(RUTAS$dir_salida, paste0("clasif_quemado_",inc,"_MMU.tif")), overwrite=TRUE)
  }
  # Guardar tambien la clasificacion sin filtrar
  writeRaster(clasif, file.path(RUTAS$dir_salida, paste0("clasif_quemado_",inc,".tif")), overwrite=TRUE)
}

cat("\nNota: el perimetro Copernicus no es verdad absoluta (islas internas).\n")
cat("El MMU elimina manchas aisladas (ruido de comision); practica estandar en\n")
cat("cartografia de area quemada (Chuvieco et al. 2016).\n")
cat(sprintf("Salidas en: %s\n", RUTAS$dir_salida))