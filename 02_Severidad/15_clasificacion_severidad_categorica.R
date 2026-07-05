# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · CLASIFICACION DE SEVERIDAD (CBI categorico)
# =============================================================================
# Espejo del analisis continuo, pero clasificando el CBI_total en 3 categorias:
#   Baja (CBI<1.0) / Moderada (1.0-2.0) / Alta (>2.0).
# Umbrales 1.0 y 2.0: tercios iguales de la escala CBI (0-3), clases
# equilibradas (80/88/67), rango compatible con literatura. NO se usan umbrales
# por incendio: el CBI es medida estandarizada y comparable (Key & Benson 2006).
#
# Estructura identica al continuo:
#   - GLOBAL: RF clasificacion + LOFOCV. Modelos: SAR, Optico, Fusion.
#   - ESTRATIFICADO por ecosistema: k-fold (principal) + LOFOCV (robustez).
#   - Metricas: accuracy, kappa (Cohen), matriz de confusion.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(dplyr)
})

RUTAS <- list(
  csv_cbi     = "G:/TFM/_DATOS/CBI_incendios.csv",
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  rds_opt_eco = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/optico_por_ecosistema_v3.rds",
  mf_raster   = "G:/TFM/_CARTOGRAFIA/FORESTAL/MF_raster.tif",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
INCENDIOS_EXCLUIR <- c("Cuevas","SanBartolome","Ciperez")
N_ARBOLES <- 500
K_FOLDS   <- 5
SEED <- 42
MF_COD <- c("1"="arbolado", "2"="matorral")
CORTES <- c(-0.01, 1.0, 2.0, 3.01)
ETIQUETAS <- c("Baja","Moderada","Alta")

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
opt_eco <- readRDS(RUTAS$rds_opt_eco)
dataset    <- rds$dataset
bandas_sar_dep <- coln$bandas_depuradas
bandas_opt_dep <- opt_eco$bandas_opt_dep
bandas_fusion  <- c(bandas_sar_dep, bandas_opt_dep)

# Variable categorica
dataset$sev_clase <- cut(dataset$CBI_total, breaks=CORTES, labels=ETIQUETAS)

cat("=============================================================\n")
cat("CLASIFICACION DE SEVERIDAD (CBI 3 categorias)\n")
cat("=============================================================\n")
cat("Distribucion de clases:\n"); print(table(dataset$sev_clase))

# --- Asignar ecosistema ---
dataset$eco <- NA_character_
cbi_csv <- read.csv(RUTAS$csv_cbi)
cbi_csv <- cbi_csv[!cbi_csv$Incendio %in% INCENDIOS_EXCLUIR, ]
idx <- match(paste(round(dataset$x_utm), round(dataset$y_utm)),
             paste(round(cbi_csv$X), round(cbi_csv$Y)))
campo <- dataset$origen == "campo"
dataset$eco[campo & !is.na(idx[campo])] <- cbi_csv$Ecosistema_L3[idx[campo][!is.na(idx[campo])]]
mf <- rast(RUTAS$mf_raster)
cero <- dataset$origen == "cero"
pts_cero <- vect(dataset[cero,], geom=c("x_utm","y_utm"), crs="EPSG:25830")
dataset$eco[cero] <- MF_COD[as.character(terra::extract(mf, pts_cero)[,2])]

# --- Metricas de clasificacion (accuracy, kappa) ---
calc_clasif <- function(obs, pred) {
  obs <- factor(obs, levels=ETIQUETAS); pred <- factor(pred, levels=ETIQUETAS)
  cm <- table(obs, pred)
  acc <- sum(diag(cm)) / sum(cm)
  # Kappa de Cohen
  n <- sum(cm); p_obs <- acc
  p_exp <- sum(rowSums(cm) * colSums(cm)) / n^2
  kappa <- (p_obs - p_exp) / (1 - p_exp)
  list(accuracy=acc, kappa=kappa, cm=cm)
}

# --- LOFOCV clasificacion ---
lofocv_clf <- function(dat, bandas, etiqueta, mostrar_cm=FALSE) {
  form <- as.formula(paste("sev_clase ~", paste(bandas, collapse=" + ")))
  pred_all <- data.frame()
  for (inc in unique(dat$Incendio)) {
    tr <- dat[dat$Incendio!=inc,]; te <- dat[dat$Incendio==inc,]
    if (nrow(te)==0) next
    set.seed(SEED)
    rf <- ranger(form, data=tr, num.trees=N_ARBOLES, seed=SEED)
    pred <- predict(rf, te)$predictions
    pred_all <- rbind(pred_all, data.frame(obs=as.character(te$sev_clase),
                                           pred=as.character(pred)))
  }
  r <- calc_clasif(pred_all$obs, pred_all$pred)
  cat(sprintf("\n=== %s (LOFOCV, n=%d) ===\n", etiqueta, nrow(pred_all)))
  cat(sprintf("  Accuracy=%.3f | Kappa=%.3f\n", r$accuracy, r$kappa))
  if (mostrar_cm) { cat("  Matriz de confusion (obs filas / pred columnas):\n"); print(r$cm) }
  data.frame(modelo=etiqueta, n=nrow(pred_all),
             accuracy=round(r$accuracy,3), kappa=round(r$kappa,3))
}

# --- k-fold clasificacion ---
kfold_clf <- function(dat, bandas, etiqueta) {
  set.seed(SEED)
  form <- as.formula(paste("sev_clase ~", paste(bandas, collapse=" + ")))
  folds <- sample(rep(1:K_FOLDS, length.out=nrow(dat)))
  pred_all <- data.frame()
  for (k in 1:K_FOLDS) {
    tr <- dat[folds!=k,]; te <- dat[folds==k,]
    rf <- ranger(form, data=tr, num.trees=N_ARBOLES, seed=SEED)
    pred <- predict(rf, te)$predictions
    pred_all <- rbind(pred_all, data.frame(obs=as.character(te$sev_clase),
                                           pred=as.character(pred)))
  }
  r <- calc_clasif(pred_all$obs, pred_all$pred)
  cat(sprintf("  %s (k-fold, n=%d): Accuracy=%.3f | Kappa=%.3f\n",
              etiqueta, nrow(dat), r$accuracy, r$kappa))
  data.frame(modelo=etiqueta, n=nrow(dat),
             accuracy=round(r$accuracy,3), kappa=round(r$kappa,3))
}

# =============================================================================
# 1. GLOBAL (LOFOCV) - 3 modelos
# =============================================================================
cat("\n=============================================================\n")
cat("1. GLOBAL (LOFOCV) · SAR vs Optico vs Fusion\n")
cat("=============================================================\n")
glob <- rbind(
  lofocv_clf(dataset, bandas_sar_dep, "SAR",    mostrar_cm=TRUE),
  lofocv_clf(dataset, bandas_opt_dep, "Optico", mostrar_cm=TRUE),
  lofocv_clf(dataset, bandas_fusion,  "Fusion", mostrar_cm=TRUE)
)
cat("\n--- Resumen global ---\n")
print(glob, row.names=FALSE)

# =============================================================================
# 2. ESTRATIFICADO POR ECOSISTEMA (k-fold principal + LOFOCV robustez)
# =============================================================================
cat("\n=============================================================\n")
cat("2. POR ECOSISTEMA · k-fold (principal) + LOFOCV (robustez)\n")
cat("=============================================================\n")
dat_eco <- dataset[!is.na(dataset$eco),]

for (e in c("arbolado","matorral")) {
  cat(sprintf("\n----- %s -----\n", toupper(e)))
  d <- dat_eco[dat_eco$eco==e,]
  cat("k-fold:\n")
  kf <- rbind(
    kfold_clf(d, bandas_sar_dep, "SAR"),
    kfold_clf(d, bandas_opt_dep, "Optico"),
    kfold_clf(d, bandas_fusion,  "Fusion")
  )
  cat("LOFOCV (robustez):\n")
  lo <- rbind(
    lofocv_clf(d, bandas_sar_dep, "SAR"),
    lofocv_clf(d, bandas_opt_dep, "Optico"),
    lofocv_clf(d, bandas_fusion,  "Fusion")
  )
  assign(paste0("kf_",e), kf); assign(paste0("lo_",e), lo)
}

cat("\n=============================================================\n")
cat("RESUMEN FINAL · KAPPA POR ECOSISTEMA\n")
cat("=============================================================\n")
cat("                    ARBOLADO            MATORRAL\n")
cat("                 k-fold  LOFOCV      k-fold  LOFOCV\n")
for (i in 1:3) {
  cat(sprintf("  %-8s        %.3f   %.3f       %.3f   %.3f\n",
              kf_arbolado$modelo[i], kf_arbolado$kappa[i], lo_arbolado$kappa[i],
              kf_matorral$kappa[i], lo_matorral$kappa[i]))
}

saveRDS(list(global=glob, kf_arbolado=kf_arbolado, lo_arbolado=lo_arbolado,
             kf_matorral=kf_matorral, lo_matorral=lo_matorral,
             cortes=CORTES, distribucion_clases=table(dataset$sev_clase)),
        file.path(RUTAS$dir_salida, "clasificacion_severidad_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "clasificacion_severidad_v3.rds")))