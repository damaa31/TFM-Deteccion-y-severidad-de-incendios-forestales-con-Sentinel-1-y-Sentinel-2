# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · ANÁLISIS POR ECOSISTEMA (arbolado vs matorral)
# =============================================================================
# Pregunta física: ¿el SAR (banda C) gradúa mejor la severidad en matorral que
# en arbolado? En arbolado denso la banda C satura en el dosel; en matorral, sin
# dosel, el cambio estructural es más directo y proporcional a la severidad.
#
#   OPCIÓN C — Univariante Spearman por ecosistema (solo parcelas de campo, 185).
#   OPCIÓN B — Modelos RF con k-fold por ecosistema (campo + puntos CBI=0
#              clasificados con la capa forestal MF_raster). Se usa k-fold (no
#              LOFOCV) por el menor tamaño muestral al estratificar.
#
# Capa forestal: MF_raster.tif (10 m, EPSG:25830). Codificacion: 1=arbolado,
# 2=matorral.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(dplyr)
})

RUTAS <- list(
  csv_cbi     = "G:/TFM/_DATOS/CBI_incendios.csv",
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  mf_raster   = "G:/TFM/_CARTOGRAFIA/FORESTAL/MF_raster.tif",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
INCENDIOS_EXCLUIR <- c("Cuevas","SanBartolome","Ciperez")
N_ARBOLES <- 500
K_FOLDS   <- 5
SEED <- 42
MF_COD <- c("1"="arbolado", "2"="matorral")

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
dataset    <- rds$dataset
bandas_sar <- rds$bandas_sar
bandas_dep <- coln$bandas_depuradas

cat("=============================================================\n")
cat("ANALISIS POR ECOSISTEMA (arbolado vs matorral)\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# Asignar ecosistema
# -----------------------------------------------------------------------------
dataset$eco <- NA_character_

cbi_csv <- read.csv(RUTAS$csv_cbi)
cbi_csv <- cbi_csv[!cbi_csv$Incendio %in% INCENDIOS_EXCLUIR, ]
key_ds  <- paste(round(dataset$x_utm), round(dataset$y_utm))
key_csv <- paste(round(cbi_csv$X), round(cbi_csv$Y))
idx <- match(key_ds, key_csv)
campo <- dataset$origen == "campo"
dataset$eco[campo & !is.na(idx[campo])] <- cbi_csv$Ecosistema_L3[idx[campo][!is.na(idx[campo])]]

mf <- rast(RUTAS$mf_raster)
cero <- dataset$origen == "cero"
pts_cero <- vect(dataset[cero, ], geom=c("x_utm","y_utm"), crs="EPSG:25830")
val_mf <- terra::extract(mf, pts_cero)[,2]
dataset$eco[cero] <- MF_COD[as.character(val_mf)]

cat(sprintf("Ecosistema asignado: %d campo, %d puntos CBI=0\n",
            sum(!is.na(dataset$eco) & campo), sum(!is.na(dataset$eco) & cero)))
cat("\n--- Distribucion por ecosistema y origen ---\n")
print(table(Ecosistema=dataset$eco, Origen=dataset$origen))

# =============================================================================
# OPCION C - UNIVARIANTE POR ECOSISTEMA (solo parcelas de campo)
# =============================================================================
cat("\n=============================================================\n")
cat("OPCION C - UNIVARIANTE SAR POR ECOSISTEMA (Spearman, solo campo)\n")
cat("=============================================================\n")
dat_campo <- dataset[campo & !is.na(dataset$eco), ]
cat(sprintf("Parcelas: %d arbolado, %d matorral\n",
            sum(dat_campo$eco=="arbolado"), sum(dat_campo$eco=="matorral")))

univ_eco <- data.frame()
for (b in bandas_sar) {
  rho_a <- cor(dat_campo[[b]][dat_campo$eco=="arbolado"],
               dat_campo$CBI_total[dat_campo$eco=="arbolado"], method="spearman")
  rho_m <- cor(dat_campo[[b]][dat_campo$eco=="matorral"],
               dat_campo$CBI_total[dat_campo$eco=="matorral"], method="spearman")
  univ_eco <- rbind(univ_eco, data.frame(
    Banda=b, rho_arbolado=round(rho_a,3), rho_matorral=round(rho_m,3),
    dif_abs=round(abs(rho_m)-abs(rho_a),3)))
}
univ_eco <- univ_eco[order(-abs(univ_eco$rho_matorral)), ]
cat("\n(dif_abs>0 = la banda correlaciona mejor con CBI en matorral)\n")
print(univ_eco, row.names = FALSE)
cat(sprintf("\n  |rho| medio arbolado: %.3f | matorral: %.3f\n",
            mean(abs(univ_eco$rho_arbolado)), mean(abs(univ_eco$rho_matorral))))
cat(sprintf("  Bandas mejores en matorral: %d de %d\n",
            sum(univ_eco$dif_abs>0), nrow(univ_eco)))

# =============================================================================
# OPCION B - MODELOS RF POR ECOSISTEMA (k-fold, campo + ceros)
# =============================================================================
cat("\n=============================================================\n")
cat("OPCION B - MODELO SAR POR ECOSISTEMA (RF + k-fold)\n")
cat("=============================================================\n")
cat(sprintf("Validacion: %d-fold (no LOFOCV) por menor tamano muestral.\n", K_FOLDS))

calc_metricas <- function(obs, pred) {
  ss_res <- sum((obs-pred)^2); ss_tot <- sum((obs-mean(obs))^2)
  data.frame(r2=1-ss_res/ss_tot, r_pearson=cor(obs,pred),
             rmse=sqrt(mean((obs-pred)^2)), mae=mean(abs(obs-pred)))
}
modelo_kfold <- function(dat, bandas, etiqueta) {
  set.seed(SEED)
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  folds <- sample(rep(1:K_FOLDS, length.out=nrow(dat)))
  pred_all <- data.frame()
  for (k in 1:K_FOLDS) {
    train <- dat[folds!=k, ]; test <- dat[folds==k, ]
    rf <- ranger(form, data=train, num.trees=N_ARBOLES, seed=SEED)
    pred <- pmin(pmax(predict(rf,test)$predictions,0),3)
    pred_all <- rbind(pred_all, data.frame(obs=test$CBI_total, pred=pred))
  }
  m <- calc_metricas(pred_all$obs, pred_all$pred)
  m$grupo <- etiqueta; m$n <- nrow(dat)
  cat(sprintf("\n=== %s (n=%d) ===\n", etiqueta, nrow(dat)))
  cat(sprintf("  R2=%.3f | r=%.3f | RMSE=%.3f | MAE=%.3f\n",
              m$r2, m$r_pearson, m$rmse, m$mae))
  m
}

dat_eco <- dataset[!is.na(dataset$eco), ]
res_b <- rbind(
  modelo_kfold(dat_eco[dat_eco$eco=="arbolado", ], bandas_dep, "SAR arbolado"),
  modelo_kfold(dat_eco[dat_eco$eco=="matorral", ], bandas_dep, "SAR matorral")
)
cat("\n--- Comparativa modelos por ecosistema ---\n")
print(res_b[, c("grupo","n","r2","r_pearson","rmse","mae")] %>%
        mutate(across(where(is.numeric), ~round(.,3))), row.names=FALSE)

saveRDS(list(univariante_eco=univ_eco, modelos_eco=res_b, k_folds=K_FOLDS,
             distribucion=table(dataset$eco, dataset$origen)),
        file.path(RUTAS$dir_salida, "severidad_por_ecosistema_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "severidad_por_ecosistema_v3.rds")))
cat("Analisis por ecosistema completado.\n")



# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · ROBUSTEZ LOFOCV POR ECOSISTEMA
# =============================================================================
# Comprobacion de robustez (NO metodo principal): repite los modelos SAR por
# ecosistema con LOFOCV (Leave-One-Fire-Out) en lugar de k-fold, para reportar
# en la memoria cuanto aguanta el hallazgo bajo validacion espacial estricta.
#
# Aviso: al estratificar, varios incendios aportan muy pocas parcelas por
# ecosistema (matorral: Canalejas=0). Por eso algunos folds seran inestables o
# daran R2 negativo; es esperable y se reporta con honestidad. El k-fold del
# script 12 es el metodo principal; esto es solo el test de robustez.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(dplyr)
})

RUTAS <- list(
  csv_cbi     = "G:/TFM/_DATOS/CBI_incendios.csv",
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  mf_raster   = "G:/TFM/_CARTOGRAFIA/FORESTAL/MF_raster.tif",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
INCENDIOS_EXCLUIR <- c("Cuevas","SanBartolome","Ciperez")
N_ARBOLES <- 500
SEED <- 42
MF_COD <- c("1"="arbolado", "2"="matorral")
MIN_TEST <- 4   # minimo de parcelas en el fold test para considerarlo fiable

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
dataset    <- rds$dataset
bandas_dep <- coln$bandas_depuradas

# --- Asignar ecosistema (igual que script 12) ---
dataset$eco <- NA_character_
cbi_csv <- read.csv(RUTAS$csv_cbi)
cbi_csv <- cbi_csv[!cbi_csv$Incendio %in% INCENDIOS_EXCLUIR, ]
key_ds  <- paste(round(dataset$x_utm), round(dataset$y_utm))
key_csv <- paste(round(cbi_csv$X), round(cbi_csv$Y))
idx <- match(key_ds, key_csv)
campo <- dataset$origen == "campo"
dataset$eco[campo & !is.na(idx[campo])] <- cbi_csv$Ecosistema_L3[idx[campo][!is.na(idx[campo])]]
mf <- rast(RUTAS$mf_raster)
cero <- dataset$origen == "cero"
pts_cero <- vect(dataset[cero, ], geom=c("x_utm","y_utm"), crs="EPSG:25830")
dataset$eco[cero] <- MF_COD[as.character(terra::extract(mf, pts_cero)[,2])]

calc_metricas <- function(obs, pred) {
  ss_res <- sum((obs-pred)^2); ss_tot <- sum((obs-mean(obs))^2)
  data.frame(r2=1-ss_res/ss_tot, r_pearson=cor(obs,pred),
             rmse=sqrt(mean((obs-pred)^2)), mae=mean(abs(obs-pred)))
}

lofocv_eco <- function(dat, bandas, etiqueta) {
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  incendios <- unique(dat$Incendio)
  pred_all <- data.frame(); por_inc <- data.frame()
  for (inc in incendios) {
    train <- dat[dat$Incendio != inc, ]; test <- dat[dat$Incendio == inc, ]
    if (nrow(test) == 0) next
    set.seed(SEED)
    rf <- ranger(form, data=train, num.trees=N_ARBOLES, seed=SEED)
    pred <- pmin(pmax(predict(rf,test)$predictions,0),3)
    pred_all <- rbind(pred_all, data.frame(obs=test$CBI_total, pred=pred))
    r_inc <- if (nrow(test) >= 2) cor(test$CBI_total, pred) else NA
    por_inc <- rbind(por_inc, data.frame(Incendio=inc, n=nrow(test),
                                         r=round(r_inc,3)))
  }
  m <- calc_metricas(pred_all$obs, pred_all$pred)
  cat(sprintf("\n=== %s (LOFOCV, n=%d) ===\n", etiqueta, nrow(pred_all)))
  cat(sprintf("  R2=%.3f | r=%.3f | RMSE=%.3f | MAE=%.3f\n",
              m$r2, m$r_pearson, m$rmse, m$mae))
  cat("  Por incendio (r de Pearson obs-pred):\n")
  print(por_inc, row.names=FALSE)
  list(global=m, por_incendio=por_inc, etiqueta=etiqueta)
}

cat("=============================================================\n")
cat("ROBUSTEZ LOFOCV POR ECOSISTEMA (test, no metodo principal)\n")
cat("=============================================================\n")

dat_eco <- dataset[!is.na(dataset$eco), ]
cat("\n--- Parcelas por incendio y ecosistema ---\n")
print(table(dat_eco$Incendio, dat_eco$eco))

r_arb <- lofocv_eco(dat_eco[dat_eco$eco=="arbolado", ], bandas_dep, "SAR arbolado")
r_mat <- lofocv_eco(dat_eco[dat_eco$eco=="matorral", ], bandas_dep, "SAR matorral")

cat("\n=============================================================\n")
cat("COMPARATIVA k-fold vs LOFOCV (para la memoria)\n")
cat("=============================================================\n")
cat("  (k-fold del script 12: arbolado R2=0.065 | matorral R2=0.423)\n")
cat(sprintf("  LOFOCV          : arbolado R2=%.3f | matorral R2=%.3f\n",
            r_arb$global$r2, r_mat$global$r2))
cat("\nLectura: si en LOFOCV el matorral sigue por encima del arbolado, el\n")
cat("hallazgo (SAR gradua severidad en matorral, no en arbolado) es robusto a\n")
cat("la validacion espacial, aunque el R2 absoluto baje por el menor tamano.\n")

saveRDS(list(lofocv_arbolado=r_arb, lofocv_matorral=r_mat),
        file.path(RUTAS$dir_salida, "ecosistema_lofocv_robustez_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "ecosistema_lofocv_robustez_v3.rds")))