# =============================================================================
# TFM SENTINEL-1 · MODELO OPTICO DE CONTROL EN DETECCION (Objetivo 1)
# =============================================================================
# Replica el modelo de deteccion SAR pero con las bandas OPTICAS (Sentinel-2),
# extraidas en los MISMOS 1.265 puntos fotointerpretados. Mismo diseno:
#   - Random Forest (ranger), 500 arboles, probabilidad.
#   - LOFOCV (deja un incendio fuera).
#   - Umbral calibrado por F1 sobre OOB del train (el test no interviene).
# Sirve como REFERENCIA OPTICA frente a la que se interpreta el SAR en deteccion.
# Salida: resultados_deteccion_OPTICO_v3.rds
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_modelado_v3.rds",  # tiene los puntos (x_utm,y_utm,quemado,Incendio)
  dir_s2      = "G:/TFM/_DATOS/EXPORT_S2_TFM",                                  # TIF opticos por incendio
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/"
)

# Bandas opticas espectrales a usar como predictores (SIN topografia, coherente con deteccion).
# Son las mismas que exporta tu script de S2; ajusta si algun nombre difiere en el TIF.
BANDAS_OPT <- c("dNBR","RBR","RdNBR","dNDMI","dBAIS2","dNBRplus","dNDVI",
                "pre_NBR","pre_NDMI","pre_BAIS2","pre_NDVI")

N_TREES <- 500
SEED    <- 42
set.seed(SEED)

cat("=============================================================\n")
cat("MODELO OPTICO DE CONTROL · DETECCION (RF + LOFOCV)\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# 1. Cargar puntos de deteccion (los mismos del modelo SAR)
# -----------------------------------------------------------------------------
rds <- readRDS(RUTAS$rds_dataset)
dataset <- rds$dataset
dataset <- dataset[dataset$Incendio != "Ciperez" & !is.na(dataset$quemado), ]
cat(sprintf("Puntos de deteccion: %d\n", nrow(dataset)))

# -----------------------------------------------------------------------------
# 2. Extraer las bandas opticas de cada TIF S2 en los puntos de su incendio
# -----------------------------------------------------------------------------
buscar_tif <- function(d, inc){
  t <- list.files(d, pattern="\\.tif$", full.names=TRUE)
  m <- t[grepl(paste0(inc,"\\.tif"), basename(t))]
  if (length(m)==0) NA_character_ else m[1]
}

incendios <- unique(dataset$Incendio)
ext_list <- list()

for (inc in incendios) {
  tif <- buscar_tif(RUTAS$dir_s2, inc)
  pts_inc <- dataset[dataset$Incendio == inc, ]
  if (is.na(tif)) { cat(sprintf("  [%s] sin TIF S2, se omite\n", inc)); next }
  r <- rast(tif)
  bandas_disp <- intersect(BANDAS_OPT, names(r))
  if (length(bandas_disp) < length(BANDAS_OPT)) {
    faltan <- setdiff(BANDAS_OPT, names(r))
    cat(sprintf("  [%s] OJO faltan bandas en TIF: %s\n", inc, paste(faltan, collapse=", ")))
  }
  r <- r[[bandas_disp]]
  v <- vect(pts_inc, geom=c("x_utm","y_utm"), crs="EPSG:25830")
  ex <- terra::extract(r, v)[, -1, drop=FALSE]   # quita columna ID
  ex$punto_id <- pts_inc$punto_id
  ex$Incendio <- pts_inc$Incendio
  ex$quemado  <- pts_inc$quemado
  ext_list[[inc]] <- ex
  cat(sprintf("  [%s] extraidos %d puntos x %d bandas opticas\n", inc, nrow(ex), length(bandas_disp)))
}

opt <- bind_rows(ext_list)
opt$quemado <- as.factor(opt$quemado)
# Quitar puntos sin dato optico (p.ej. fuera de cobertura / nube)
antes <- nrow(opt)
opt <- opt[complete.cases(opt[, BANDAS_OPT]), ]
cat(sprintf("\nPuntos con optico completo: %d (descartados %d por NA)\n", nrow(opt), antes-nrow(opt)))

# -----------------------------------------------------------------------------
# 3. Funciones (identicas al modelo SAR para comparabilidad)
# -----------------------------------------------------------------------------
calc_f1 <- function(obs, pred){
  obs<-as.integer(as.character(obs)); pred<-as.integer(as.character(pred))
  tp<-sum(obs==1&pred==1); fp<-sum(obs==0&pred==1); fn<-sum(obs==1&pred==0)
  if(tp+fp==0||tp+fn==0) return(0)
  pr<-tp/(tp+fp); rc<-tp/(tp+fn); if(pr+rc==0) return(0); 2*pr*rc/(pr+rc)
}
calibrar_umbral <- function(probs, obs, grid=seq(0.10,0.90,0.01)){
  obs<-as.integer(as.character(obs))
  f1s<-sapply(grid, function(th) calc_f1(obs, as.integer(probs>=th)))
  list(umbral=grid[which.max(f1s)], f1=max(f1s))
}
metricas <- function(obs, pred){
  obs<-as.integer(as.character(obs)); pred<-as.integer(as.character(pred))
  tp<-sum(obs==1&pred==1); fp<-sum(obs==0&pred==1); fn<-sum(obs==1&pred==0); tn<-sum(obs==0&pred==0)
  n<-tp+fp+fn+tn; acc<-(tp+tn)/n
  pr<-ifelse(tp+fp==0,NA,tp/(tp+fp)); rc<-ifelse(tp+fn==0,NA,tp/(tp+fn))
  sp<-ifelse(tn+fp==0,NA,tn/(tn+fp))
  f1<-ifelse(is.na(pr)||is.na(rc)||(pr+rc)==0,NA,2*pr*rc/(pr+rc))
  iou<-ifelse(tp+fp+fn==0,NA,tp/(tp+fp+fn))
  om<-ifelse(tp+fn==0,NA,fn/(tp+fn)); com<-ifelse(tp+fp==0,NA,fp/(tp+fp))
  pe<-((tp+fp)*(tp+fn)+(fn+tn)*(fp+tn))/(n^2); kappa<-(acc-pe)/(1-pe)
  data.frame(n=n, accuracy=acc, dice_f1=f1, iou=iou, precision=pr, recall=rc,
             especificidad=sp, omision=om, comision=com, kappa=kappa)
}

# -----------------------------------------------------------------------------
# 4. LOFOCV con las bandas opticas
# -----------------------------------------------------------------------------
metricas_incendio <- data.frame(); predicciones <- data.frame(); umbrales <- data.frame()

for (id in unique(opt$Incendio)) {
  train <- opt[opt$Incendio != id, ]
  test  <- opt[opt$Incendio == id, ]
  rf <- ranger(as.formula(paste("quemado ~", paste(BANDAS_OPT, collapse=" + "))),
               data=train, num.trees=N_TREES, probability=TRUE, seed=SEED)
  cal <- calibrar_umbral(rf$predictions[,"1"], train$quemado)
  probs_test <- predict(rf, data=test)$predictions[,"1"]
  pred_test  <- as.integer(probs_test >= cal$umbral)
  m <- metricas(test$quemado, pred_test); m$Incendio<-id; m$umbral<-cal$umbral
  metricas_incendio <- bind_rows(metricas_incendio, m)
  umbrales <- bind_rows(umbrales, data.frame(Incendio=id, umbral=cal$umbral, f1_oob=cal$f1))
  predicciones <- bind_rows(predicciones, data.frame(
    Incendio=id, punto_id=test$punto_id, obs=as.integer(as.character(test$quemado)),
    prob=probs_test, pred=pred_test))
  cat(sprintf("  [%s] Dice=%.3f Recall=%.3f Comision=%.3f Kappa=%.3f (umbral %.2f)\n",
              id, m$dice_f1, m$recall, m$comision, m$kappa, cal$umbral))
}

m_global <- metricas(predicciones$obs, predicciones$pred)

cat("\n=== METRICAS GLOBALES OPTICO (deteccion) ===\n"); print(round(m_global,3))
cat("\n=== POR INCENDIO ===\n")
print(as.data.frame(metricas_incendio %>% mutate(across(where(is.numeric),~round(.,3)))), row.names=FALSE)

saveRDS(list(metricas_global=m_global, metricas_incendio=metricas_incendio,
             umbrales=umbrales, predicciones=predicciones, bandas_opt=BANDAS_OPT,
             n_trees=N_TREES, seed=SEED),
        file.path(RUTAS$dir_salida, "resultados_deteccion_OPTICO_v3.rds"))
cat(sprintf("\nGuardado: %s\n", file.path(RUTAS$dir_salida,"resultados_deteccion_OPTICO_v3.rds")))