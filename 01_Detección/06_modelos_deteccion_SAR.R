# =============================================================================
# TFM SENTINEL-1 · MÓDULO 3 (v3): MODELADO DE DETECCIÓN
# =============================================================================
# Objetivo 1: detectar zonas quemadas vs no quemadas únicamente a partir de
# las 14 bandas Sentinel-1, usando los 1.415 puntos fotointerpretados como
# verdad terreno.
#
# Diseño:
#   - Algoritmo: Random Forest (ranger), 500 árboles, mtry por defecto.
#   - Validación: Leave-One-Fire-Out Cross-Validation (LOFOCV). 9 folds.
#   - Calibración del umbral: por F1 óptimo sobre predicciones OOB del
#     training de cada fold. El incendio de test no interviene.
#   - Desbalance: se respeta la distribución natural (43/57).
#   - Predictores: las 14 bandas SAR estándar.
#
# Métricas reportadas: accuracy, Dice (F1), IoU, precisión, recall,
# especificidad, omisión, comisión, kappa de Cohen.
#
# Salida: resultados_deteccion_v3.rds
# =============================================================================

suppressPackageStartupMessages({
  library(ranger); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_modelado_v3.rds",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/"
)

N_TREES <- 500
SEED    <- 42
set.seed(SEED)

cat("=============================================================\n")
cat("MÓDULO 3 (v3): MODELADO DE DETECCIÓN (RF + LOFOCV)\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# Cargar dataset
# -----------------------------------------------------------------------------
rds <- readRDS(RUTAS$rds_dataset)
dataset    <- rds$dataset
bandas_sar <- rds$bandas_sar
dataset$quemado <- as.factor(dataset$quemado)

cat(sprintf("Dataset: %d puntos | %d quemado | %d no quemado\n",
            nrow(dataset),
            sum(dataset$quemado == "1"),
            sum(dataset$quemado == "0")))
cat(sprintf("Predictores: %d bandas SAR\n", length(bandas_sar)))
cat(sprintf("Incendios: %s\n\n", paste(unique(dataset$Incendio), collapse=", ")))

# -----------------------------------------------------------------------------
# Funciones auxiliares
# -----------------------------------------------------------------------------
calc_f1 <- function(obs, pred) {
  obs  <- as.integer(as.character(obs))
  pred <- as.integer(as.character(pred))
  tp <- sum(obs == 1 & pred == 1)
  fp <- sum(obs == 0 & pred == 1)
  fn <- sum(obs == 1 & pred == 0)
  if (tp + fp == 0 || tp + fn == 0) return(0)
  prec <- tp / (tp + fp); rec <- tp / (tp + fn)
  if (prec + rec == 0) return(0)
  2 * prec * rec / (prec + rec)
}

calibrar_umbral <- function(probs_oob, obs_oob, grid = seq(0.10, 0.90, 0.01)) {
  obs <- as.integer(as.character(obs_oob))
  f1s <- sapply(grid, function(th) calc_f1(obs, as.integer(probs_oob >= th)))
  list(umbral = grid[which.max(f1s)], f1 = max(f1s))
}

metricas <- function(obs, pred) {
  obs  <- as.integer(as.character(obs))
  pred <- as.integer(as.character(pred))
  tp <- sum(obs == 1 & pred == 1); fp <- sum(obs == 0 & pred == 1)
  fn <- sum(obs == 1 & pred == 0); tn <- sum(obs == 0 & pred == 0)
  n  <- tp + fp + fn + tn
  
  acc  <- (tp + tn) / n
  prec <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  rec  <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  spec <- ifelse(tn + fp == 0, NA, tn / (tn + fp))
  f1   <- ifelse(is.na(prec)||is.na(rec)||(prec+rec)==0, NA, 2*prec*rec/(prec+rec))
  iou  <- ifelse(tp + fp + fn == 0, NA, tp / (tp + fp + fn))
  om   <- ifelse(tp + fn == 0, NA, fn / (tp + fn))
  com  <- ifelse(tp + fp == 0, NA, fp / (tp + fp))
  
  po <- acc
  pe <- ((tp+fp)*(tp+fn) + (fn+tn)*(fp+tn)) / (n^2)
  kappa <- (po - pe) / (1 - pe)
  
  data.frame(n=n, accuracy=acc, dice_f1=f1, iou=iou,
             precision=prec, recall=rec, especificidad=spec,
             omision=om, comision=com, kappa=kappa)
}

# -----------------------------------------------------------------------------
# LOFOCV
# -----------------------------------------------------------------------------
incendios <- unique(dataset$Incendio)
predicciones_globales <- data.frame()
metricas_por_incendio <- data.frame()
umbrales_calibrados   <- data.frame()

for (id in incendios) {
  cat(sprintf("--- Fold: test = %s ---\n", id))
  
  train <- dataset[dataset$Incendio != id, ]
  test  <- dataset[dataset$Incendio == id, ]
  
  rf <- ranger(
    formula     = as.formula(paste("quemado ~", paste(bandas_sar, collapse=" + "))),
    data        = train,
    num.trees   = N_TREES,
    probability = TRUE,
    seed        = SEED
  )
  
  probs_oob_train <- rf$predictions[, "1"]
  obs_train       <- train$quemado
  cal <- calibrar_umbral(probs_oob_train, obs_train)
  cat(sprintf("  Umbral calibrado: %.2f (F1 OOB train: %.3f)\n",
              cal$umbral, cal$f1))
  
  probs_test <- predict(rf, data = test)$predictions[, "1"]
  pred_test  <- as.integer(probs_test >= cal$umbral)
  
  m <- metricas(test$quemado, pred_test)
  m$Incendio <- id
  m$umbral   <- cal$umbral
  metricas_por_incendio <- bind_rows(metricas_por_incendio, m)
  umbrales_calibrados <- bind_rows(umbrales_calibrados,
                                   data.frame(Incendio=id, umbral=cal$umbral, f1_oob_train=cal$f1))
  
  predicciones_globales <- bind_rows(predicciones_globales,
                                     data.frame(Incendio=id, punto_id=test$punto_id,
                                                obs=as.integer(as.character(test$quemado)),
                                                prob=probs_test, pred=pred_test, umbral=cal$umbral))
  
  cat(sprintf("  Test: Dice=%.3f | Recall=%.3f | Comisión=%.3f | Kappa=%.3f\n\n",
              m$dice_f1, m$recall, m$comision, m$kappa))
}

# -----------------------------------------------------------------------------
# Métricas globales
# -----------------------------------------------------------------------------
m_global <- metricas(predicciones_globales$obs, predicciones_globales$pred)

cat("=============================================================\n")
cat("RESULTADOS\n")
cat("=============================================================\n")

cat("\n--- Métricas globales (pool de todos los incendios) ---\n")
print(round(m_global, 3))

cat("\n--- Métricas por incendio ---\n")
tabla_incendios <- metricas_por_incendio %>%
  select(Incendio, n, accuracy, dice_f1, iou, precision, recall,
         especificidad, omision, comision, kappa, umbral)
print(as.data.frame(tabla_incendios %>%
                      mutate(across(where(is.numeric), ~round(.,3)))),
      row.names = FALSE)

cat("\n--- Umbrales calibrados por fold ---\n")
print(as.data.frame(umbrales_calibrados %>%
                      mutate(across(where(is.numeric), ~round(.,3)))),
      row.names = FALSE)
cat(sprintf("\nUmbral medio: %.3f | rango: [%.3f, %.3f]\n",
            mean(umbrales_calibrados$umbral),
            min(umbrales_calibrados$umbral),
            max(umbrales_calibrados$umbral)))

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
saveRDS(list(
  metricas_global    = m_global,
  metricas_incendio  = metricas_por_incendio,
  umbrales           = umbrales_calibrados,
  predicciones       = predicciones_globales,
  bandas_sar         = bandas_sar,
  n_trees            = N_TREES,
  seed               = SEED
), file.path(RUTAS$dir_salida, "resultados_deteccion_v3.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "resultados_deteccion_v3.rds")))
cat("Módulo 3 v3 completado.\n")



# ----------------------------------------------------------------------------------------------------------------------------------------------





# =============================================================================
# TFM SENTINEL-1 · DIAGNÓSTICO DE INCENDIOS DÉBILES (Fasgar, Barniedo)
# =============================================================================
# Inspecciona los dos incendios con menor kappa del modelo SAR final:
#   - Matrices de confusión detalladas.
#   - Distribución de las bandas SAR más informativas (top Boruta).
#   - Comparación con los incendios "fuertes" (Porto, Llamas).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

RUTAS <- list(
  rds_dataset     = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_modelado_v3.rds",
  rds_resultados  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/resultados_deteccion_v3.rds"
)

INCENDIOS_DEBILES <- c("Fasgar", "Barniedo")
INCENDIOS_FUERTES <- c("Porto", "Llamas")

cat("=============================================================\n")
cat("DIAGNÓSTICO DE INCENDIOS DÉBILES\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# Cargar
# -----------------------------------------------------------------------------
ds  <- readRDS(RUTAS$rds_dataset)
res <- readRDS(RUTAS$rds_resultados)
dataset       <- ds$dataset
bandas_sar    <- ds$bandas_sar
predicciones  <- res$predicciones

# Bandas más informativas según Boruta (top 4) + textura más débil
bandas_top <- c("Delta_RFDI","Delta_RVI","RBR_VH_Lin","Delta_VH_dB","STD70_DeltaVH")

# -----------------------------------------------------------------------------
# Matrices de confusión detalladas
# -----------------------------------------------------------------------------
mostrar_matriz <- function(id) {
  pp <- predicciones[predicciones$Incendio == id, ]
  if (nrow(pp) == 0) { cat(sprintf("Sin predicciones para %s\n", id)); return() }
  
  tab <- table(Obs = pp$obs, Pred = pp$pred)
  cat(sprintf("\n--- %s (n = %d, umbral = %.2f) ---\n",
              id, nrow(pp), unique(pp$umbral)))
  print(tab)
  
  # Falsos positivos y falsos negativos en términos relativos
  tp <- sum(pp$obs == 1 & pp$pred == 1)
  fp <- sum(pp$obs == 0 & pp$pred == 1)
  fn <- sum(pp$obs == 1 & pp$pred == 0)
  tn <- sum(pp$obs == 0 & pp$pred == 0)
  cat(sprintf("  TP=%d | FP=%d | FN=%d | TN=%d\n", tp, fp, fn, tn))
  cat(sprintf("  Falsos positivos: %.1f%% de las predicciones 'quemado'\n",
              100*fp/max(tp+fp,1)))
  cat(sprintf("  Falsos negativos: %.1f%% de los quemados reales\n",
              100*fn/max(tp+fn,1)))
  
  # Distribución de probabilidades
  cat(sprintf("  Probabilidades: media obs=1 = %.3f | media obs=0 = %.3f\n",
              mean(pp$prob[pp$obs==1]), mean(pp$prob[pp$obs==0])))
}

cat("\n############### INCENDIOS DÉBILES ###############\n")
for (id in INCENDIOS_DEBILES) mostrar_matriz(id)

cat("\n############### INCENDIOS FUERTES (referencia) ###############\n")
for (id in INCENDIOS_FUERTES) mostrar_matriz(id)

# -----------------------------------------------------------------------------
# Rangos de las bandas SAR top por incendio
# -----------------------------------------------------------------------------
cat("\n\n=============================================================\n")
cat("RANGOS DE LAS BANDAS SAR CLAVE\n")
cat("=============================================================\n")
cat("Comparativa débiles (Fasgar, Barniedo) vs fuertes (Porto, Llamas)\n\n")

incendios_comp <- c(INCENDIOS_DEBILES, INCENDIOS_FUERTES)

for (banda in bandas_top) {
  cat(sprintf("\n--- %s ---\n", banda))
  resumen <- dataset %>%
    filter(Incendio %in% incendios_comp) %>%
    group_by(Incendio, quemado) %>%
    summarise(
      media = round(mean(.data[[banda]], na.rm = TRUE), 3),
      sd    = round(sd(.data[[banda]],   na.rm = TRUE), 3),
      mediana = round(median(.data[[banda]], na.rm = TRUE), 3),
      q25 = round(quantile(.data[[banda]], 0.25, na.rm = TRUE), 3),
      q75 = round(quantile(.data[[banda]], 0.75, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    mutate(grupo = ifelse(Incendio %in% INCENDIOS_DEBILES, "débil", "fuerte"))
  print(as.data.frame(resumen), row.names = FALSE)
  
  # Separabilidad: diferencia entre medias quemado vs no quemado por incendio
  sep <- dataset %>%
    filter(Incendio %in% incendios_comp) %>%
    group_by(Incendio) %>%
    summarise(
      mu_quemado    = mean(.data[[banda]][quemado == 1], na.rm = TRUE),
      mu_no_quemado = mean(.data[[banda]][quemado == 0], na.rm = TRUE),
      diferencia    = round(mu_quemado - mu_no_quemado, 3),
      sd_pool       = sd(.data[[banda]], na.rm = TRUE),
      d_cohen       = round((mu_quemado - mu_no_quemado) / sd_pool, 3),
      .groups = "drop"
    )
  cat("\nSeparabilidad quemado vs no quemado (d de Cohen):\n")
  print(as.data.frame(sep %>% select(Incendio, diferencia, d_cohen)),
        row.names = FALSE)
}

cat("\n\nNota: |d de Cohen| > 0.8 = separabilidad fuerte; 0.5-0.8 = media; <0.5 = baja\n")
cat("Si los incendios débiles tienen d de Cohen muy bajos en bandas clave,\n")
cat("eso explica por qué el modelo falla en ellos: la señal SAR no separa bien las clases.\n")

