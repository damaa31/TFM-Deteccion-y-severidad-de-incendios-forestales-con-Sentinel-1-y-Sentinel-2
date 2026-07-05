# =============================================================================
# TFM SENTINEL-1 · BORUTA: SELECCIÓN DE VARIABLES SAR
# =============================================================================
# Aplica el algoritmo Boruta sobre el dataset SAR de detección para identificar
# qué bandas son confirmadas como informativas, tentativas o rechazadas.
#
# Referencia: Kursa, M. B., & Rudnicki, W. R. (2010). Feature selection with
# the Boruta package. Journal of Statistical Software, 36(11), 1-13.
#
# Lógica:
#   - Boruta compara la importancia de cada variable real con la de versiones
#     aleatorizadas ("shadow features"). Sólo confirma como informativas
#     aquellas que superan significativamente a las versiones aleatorias.
#   - Robusto al ruido y a la multicolinealidad, recomendado en estudios de
#     teledetección con stacks de muchas variables.
#
# Salida: resultados_boruta_v3.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Boruta); library(ranger); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_modelado_v3.rds",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/"
)

SEED <- 42
set.seed(SEED)

cat("=============================================================\n")
cat("BORUTA · SELECCIÓN DE VARIABLES SAR\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# Cargar dataset
# -----------------------------------------------------------------------------
rds <- readRDS(RUTAS$rds_dataset)
dataset    <- rds$dataset
bandas_sar <- rds$bandas_sar
dataset$quemado <- as.factor(dataset$quemado)

cat(sprintf("Dataset cargado: %d puntos | %d predictores\n",
            nrow(dataset), length(bandas_sar)))
cat(sprintf("Incendios: %s\n\n", paste(unique(dataset$Incendio), collapse=", ")))

# -----------------------------------------------------------------------------
# Boruta sobre el dataset completo
# -----------------------------------------------------------------------------
cat("Ejecutando Boruta (puede tardar varios minutos)...\n")
formula_boruta <- as.formula(paste("quemado ~", paste(bandas_sar, collapse=" + ")))

boruta_result <- Boruta(
  formula     = formula_boruta,
  data        = dataset,
  maxRuns     = 100,
  doTrace     = 2,
  getImp      = getImpRfZ   # importancia basada en ranger (rápido)
)

cat("\n=============================================================\n")
cat("RESULTADO DE BORUTA\n")
cat("=============================================================\n")
print(boruta_result)

# Resumen estructurado por decisión
decision <- boruta_result$finalDecision
tabla_decision <- data.frame(
  Variable = names(decision),
  Decisión = as.character(decision),
  Imp_med  = round(attStats(boruta_result)$medianImp, 3),
  Imp_min  = round(attStats(boruta_result)$minImp, 3),
  Imp_max  = round(attStats(boruta_result)$maxImp, 3)
)
tabla_decision <- tabla_decision[order(-tabla_decision$Imp_med), ]
cat("\n--- Detalle por variable (ordenado por importancia media) ---\n")
print(tabla_decision, row.names = FALSE)

# Listas explícitas
confirmadas <- names(decision)[decision == "Confirmed"]
tentativas  <- names(decision)[decision == "Tentative"]
rechazadas  <- names(decision)[decision == "Rejected"]

cat(sprintf("\nConfirmadas (%d): %s\n",
            length(confirmadas), paste(confirmadas, collapse = ", ")))
cat(sprintf("Tentativas (%d): %s\n",
            length(tentativas), paste(tentativas, collapse = ", ")))
cat(sprintf("Rechazadas (%d): %s\n",
            length(rechazadas), paste(rechazadas, collapse = ", ")))

# -----------------------------------------------------------------------------
# Resolver las tentativas (TentativeRoughFix asigna decisión por proximidad)
# -----------------------------------------------------------------------------
if (length(tentativas) > 0) {
  cat("\nResolviendo variables tentativas con TentativeRoughFix...\n")
  boruta_fixed <- TentativeRoughFix(boruta_result)
  decision_fix <- boruta_fixed$finalDecision
  confirmadas_fin <- names(decision_fix)[decision_fix == "Confirmed"]
  cat(sprintf("Confirmadas tras resolver tentativas (%d): %s\n",
              length(confirmadas_fin), paste(confirmadas_fin, collapse = ", ")))
} else {
  confirmadas_fin <- confirmadas
}

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
saveRDS(list(
  boruta            = boruta_result,
  tabla_decision    = tabla_decision,
  confirmadas       = confirmadas,
  tentativas        = tentativas,
  rechazadas        = rechazadas,
  confirmadas_final = confirmadas_fin,
  seed              = SEED
), file.path(RUTAS$dir_salida, "resultados_boruta_v3.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "resultados_boruta_v3.rds")))
cat("Selección de variables Boruta completada.\n")