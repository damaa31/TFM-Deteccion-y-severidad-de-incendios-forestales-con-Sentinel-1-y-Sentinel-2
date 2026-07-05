# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · MODELO SAR DE SEVERIDAD (núcleo)
# =============================================================================
# Random Forest de REGRESIÓN para predecir CBI_total a partir de las bandas SAR
# depuradas por colinealidad. Validación LOFOCV (Leave-One-Fire-Out), coherente
# con el Objetivo 1. Métricas: R2, RMSE, MAE (global y por incendio).
#
# Literatura: Tanase et al. (2014) -> RF para severidad SAR; Roberts et al.
# (2017) -> validación espacial (LOFOCV) para evitar optimismo por
# autocorrelación.
# =============================================================================

suppressPackageStartupMessages({
  library(ranger); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
N_ARBOLES <- 500
SEED <- 42
set.seed(SEED)

cat("=============================================================\n")
cat("OBJETIVO 2 · MODELO SAR DE SEVERIDAD (RF + LOFOCV)\n")
cat("=============================================================\n")

# Cargar dataset y bandas depuradas
rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
dataset    <- rds$dataset
bandas_sar <- coln$bandas_depuradas

cat(sprintf("Dataset: %d registros | %d bandas SAR depuradas\n",
            nrow(dataset), length(bandas_sar)))
cat("Bandas:", paste(bandas_sar, collapse=", "), "\n")

incendios <- unique(dataset$Incendio)

# -----------------------------------------------------------------------------
# Métricas de regresión
# -----------------------------------------------------------------------------
calc_metricas <- function(obs, pred) {
  rmse <- sqrt(mean((obs - pred)^2))
  mae  <- mean(abs(obs - pred))
  # R2 como 1 - SSres/SStot (coeficiente de determinación)
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - ss_res/ss_tot
  # r de Pearson obs-pred (correlación)
  r_pearson <- cor(obs, pred)
  data.frame(r2 = r2, r_pearson = r_pearson, rmse = rmse, mae = mae)
}

# -----------------------------------------------------------------------------
# LOFOCV
# -----------------------------------------------------------------------------
cat("\n--- LOFOCV (Leave-One-Fire-Out) ---\n")
form <- as.formula(paste("CBI_total ~", paste(bandas_sar, collapse = " + ")))

pred_global <- data.frame()
metricas_inc <- data.frame()

for (inc_test in incendios) {
  train <- dataset[dataset$Incendio != inc_test, ]
  test  <- dataset[dataset$Incendio == inc_test, ]
  
  rf <- ranger(form, data = train, num.trees = N_ARBOLES,
               importance = "none", seed = SEED)
  pred <- predict(rf, test)$predictions
  # Recortar predicciones al rango válido de CBI [0, 3]
  pred <- pmin(pmax(pred, 0), 3)
  
  m <- calc_metricas(test$CBI_total, pred)
  m$Incendio <- inc_test
  m$n <- nrow(test)
  metricas_inc <- rbind(metricas_inc, m)
  
  pred_global <- rbind(pred_global, data.frame(
    Incendio = inc_test, obs = test$CBI_total, pred = pred))
}

# -----------------------------------------------------------------------------
# Métricas globales (pool de todas las predicciones LOFOCV)
# -----------------------------------------------------------------------------
m_global <- calc_metricas(pred_global$obs, pred_global$pred)
cat("\n--- Métricas globales (pool LOFOCV) ---\n")
cat(sprintf("  n        : %d\n", nrow(pred_global)))
cat(sprintf("  R2       : %.3f\n", m_global$r2))
cat(sprintf("  r Pearson: %.3f\n", m_global$r_pearson))
cat(sprintf("  RMSE     : %.3f (en unidades de CBI, rango 0-3)\n", m_global$rmse))
cat(sprintf("  MAE      : %.3f\n", m_global$mae))

cat("\n--- Métricas por incendio ---\n")
ti <- metricas_inc %>%
  select(Incendio, n, r2, r_pearson, rmse, mae) %>%
  mutate(across(where(is.numeric), ~round(., 3)))
print(as.data.frame(ti), row.names = FALSE)

# -----------------------------------------------------------------------------
# Modelo final sobre todo el dataset (para importancia de variables)
# -----------------------------------------------------------------------------
rf_full <- ranger(form, data = dataset, num.trees = N_ARBOLES,
                  importance = "permutation", seed = SEED)
imp <- sort(rf_full$variable.importance, decreasing = TRUE)
cat("\n--- Importancia de variables (modelo completo) ---\n")
for (nm in names(imp)) cat(sprintf("  %-16s %.4f\n", nm, imp[nm]))

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
saveRDS(list(
  metricas_global = m_global,
  metricas_por_incendio = metricas_inc,
  predicciones = pred_global,
  importancia = imp,
  bandas = bandas_sar
), file.path(RUTAS$dir_salida, "modelo_sar_severidad_v3.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "modelo_sar_severidad_v3.rds")))
cat("Modelo SAR de severidad completado.\n")