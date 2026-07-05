# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · UNIVARIANTE + COLINEALIDAD (SAR)
# =============================================================================
# 1. Correlación de Spearman de cada banda SAR con CBI_total (univariante):
#    qué índices radar capturan severidad por sí solos.
# 2. Matriz de colinealidad Spearman entre bandas SAR: identificar redundancias
#    (|rho| > 0.90) y depurar conservando, de cada par redundante, la de mayor
#    correlación con CBI_total (criterio objetivo).
#
# Literatura: Lasaponara et al. (2019), Imperatore et al. (2017) -> univariante;
#             Dormann et al. (2013) -> colinealidad.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
UMBRAL_COLIN <- 0.90

cat("=============================================================\n")
cat("OBJETIVO 2 · UNIVARIANTE + COLINEALIDAD (SAR)\n")
cat("=============================================================\n")

rds <- readRDS(RUTAS$rds_dataset)
dataset    <- rds$dataset
bandas_sar <- rds$bandas_sar

cat(sprintf("Dataset: %d registros | %d bandas SAR\n",
            nrow(dataset), length(bandas_sar)))

# -----------------------------------------------------------------------------
# 1. UNIVARIANTE: Spearman de cada banda SAR con CBI_total
# -----------------------------------------------------------------------------
cat("\n--- 1. Correlación univariante (Spearman) banda SAR vs CBI_total ---\n")
univ <- data.frame()
for (b in bandas_sar) {
  ct <- cor.test(dataset[[b]], dataset$CBI_total, method = "spearman",
                 exact = FALSE)
  univ <- rbind(univ, data.frame(
    Banda = b,
    rho   = round(as.numeric(ct$estimate), 3),
    p_val = signif(ct$p.value, 3),
    abs_rho = abs(round(as.numeric(ct$estimate), 3))
  ))
}
univ <- univ[order(-univ$abs_rho), ]
print(univ[, c("Banda","rho","p_val")], row.names = FALSE)

cat("\nLas bandas con |rho| más alto capturan mejor el gradiente de severidad\n")
cat("por sí solas. Las de p_val > 0.05 no muestran relación significativa.\n")

# -----------------------------------------------------------------------------
# 2. COLINEALIDAD: matriz Spearman entre bandas SAR
# -----------------------------------------------------------------------------
cat("\n--- 2. Matriz de colinealidad Spearman entre bandas SAR ---\n")
mat <- cor(dataset[, bandas_sar], method = "spearman", use = "complete.obs")
print(round(mat, 2))

# Identificar pares redundantes (|rho| > umbral)
cat(sprintf("\n--- Pares con |rho| > %.2f (redundantes) ---\n", UMBRAL_COLIN))
pares_red <- data.frame()
for (i in 1:(length(bandas_sar)-1)) {
  for (j in (i+1):length(bandas_sar)) {
    r <- mat[i, j]
    if (abs(r) > UMBRAL_COLIN) {
      pares_red <- rbind(pares_red, data.frame(
        Banda_A = bandas_sar[i], Banda_B = bandas_sar[j],
        rho = round(r, 3)))
    }
  }
}
if (nrow(pares_red) == 0) {
  cat("No hay pares por encima del umbral. No se elimina ninguna banda.\n")
  bandas_depuradas <- bandas_sar
} else {
  print(pares_red, row.names = FALSE)
  
  # Depuración: de cada par redundante, eliminar la de MENOR |rho| con CBI
  cat("\n--- Depuración (conservar la de mayor |rho| con CBI) ---\n")
  abs_rho_cbi <- setNames(univ$abs_rho, univ$Banda)
  eliminar <- c()
  for (k in 1:nrow(pares_red)) {
    a <- pares_red$Banda_A[k]; b <- pares_red$Banda_B[k]
    if (a %in% eliminar || b %in% eliminar) next
    # eliminar la de menor correlación con CBI
    if (abs_rho_cbi[a] >= abs_rho_cbi[b]) {
      eliminar <- c(eliminar, b)
      cat(sprintf("  %s vs %s (rho=%.2f) -> elimina %s (menor rho con CBI)\n",
                  a, b, pares_red$rho[k], b))
    } else {
      eliminar <- c(eliminar, a)
      cat(sprintf("  %s vs %s (rho=%.2f) -> elimina %s (menor rho con CBI)\n",
                  a, b, pares_red$rho[k], a))
    }
  }
  bandas_depuradas <- setdiff(bandas_sar, eliminar)
}

cat(sprintf("\n--- Bandas SAR depuradas (%d de %d) ---\n",
            length(bandas_depuradas), length(bandas_sar)))
cat(paste(bandas_depuradas, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
saveRDS(list(
  univariante = univ,
  matriz_colin = mat,
  pares_redundantes = if (exists("pares_red")) pares_red else NULL,
  bandas_depuradas = bandas_depuradas,
  umbral = UMBRAL_COLIN
), file.path(RUTAS$dir_salida, "univariante_colinealidad_sar.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "univariante_colinealidad_sar.rds")))
cat("Análisis univariante + colinealidad completado.\n")

