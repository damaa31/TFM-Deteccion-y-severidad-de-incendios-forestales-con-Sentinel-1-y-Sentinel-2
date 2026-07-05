# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · ÓPTICO  + FUSIÓN SAR+ÓPTICO
# =============================================================================
# Cierra la comparación del Objetivo 2 (diseño Opción C):
#   1. Depura las 11 bandas ópticas por colinealidad (Spearman |rho|>0.90),
#      criterio IDÉNTICO al aplicado al SAR (conservar la de mayor |rho| con CBI).
#   2. Modelo ÓPTICO (benchmark): RF + LOFOCV con ópticas depuradas.
#   3. Modelo FUSIÓN: RF + LOFOCV con SAR depurado + óptico depurado.
#   4. Tabla comparativa final: SAR vs Óptico vs Fusión.
#
# Configuración fija (idéntica al SAR): dataset con CBI=0 (235 reg), LOFOCV.
# Literatura: Belenguer-Plomer et al. (2021) -> fusión SAR+óptico.
# =============================================================================

suppressPackageStartupMessages({
  library(ranger); library(dplyr)
})

RUTAS <- list(
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  rds_colin   = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/univariante_colinealidad_sar.rds",
  rds_sar     = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/modelo_sar_severidad_v3.rds",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
N_ARBOLES <- 500
SEED <- 42
UMBRAL_COLIN <- 0.90

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
dataset    <- rds$dataset
bandas_opt <- rds$bandas_opt
bandas_sar_dep <- coln$bandas_depuradas

cat("=============================================================\n")
cat("OBJETIVO 2 · ÓPTICO + FUSIÓN\n")
cat("=============================================================\n")

calc_metricas <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2); ss_tot <- sum((obs - mean(obs))^2)
  data.frame(r2 = 1 - ss_res/ss_tot, r_pearson = cor(obs, pred),
             rmse = sqrt(mean((obs - pred)^2)), mae = mean(abs(obs - pred)))
}

correr_lofocv <- function(dat, bandas, etiqueta) {
  set.seed(SEED)
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse = " + ")))
  pred_global <- data.frame(); met_inc <- data.frame()
  for (inc_test in unique(dat$Incendio)) {
    train <- dat[dat$Incendio != inc_test, ]; test <- dat[dat$Incendio == inc_test, ]
    rf <- ranger(form, data = train, num.trees = N_ARBOLES, seed = SEED)
    pred <- pmin(pmax(predict(rf, test)$predictions, 0), 3)
    m <- calc_metricas(test$CBI_total, pred); m$Incendio <- inc_test; m$n <- nrow(test)
    met_inc <- rbind(met_inc, m)
    pred_global <- rbind(pred_global, data.frame(Incendio=inc_test,
                                                 obs=test$CBI_total, pred=pred))
  }
  mg <- calc_metricas(pred_global$obs, pred_global$pred)
  cat(sprintf("\n=== %s ===\n", etiqueta))
  cat(sprintf("  R2=%.3f | r=%.3f | RMSE=%.3f | MAE=%.3f\n",
              mg$r2, mg$r_pearson, mg$rmse, mg$mae))
  list(global=mg, por_incendio=met_inc, pred=pred_global, etiqueta=etiqueta)
}

# -----------------------------------------------------------------------------
# 1. Depuración de colinealidad del óptico (criterio idéntico al SAR)
# -----------------------------------------------------------------------------
cat("\n--- 1. Depuración colinealidad ÓPTICO (Spearman |rho|>0.90) ---\n")
# Univariante óptico (para el criterio de conservación)
univ_opt <- sapply(bandas_opt, function(b)
  abs(cor(dataset[[b]], dataset$CBI_total, method="spearman")))
mat_opt <- cor(dataset[, bandas_opt], method = "spearman", use="complete.obs")

pares <- data.frame()
for (i in 1:(length(bandas_opt)-1)) for (j in (i+1):length(bandas_opt)) {
  if (abs(mat_opt[i,j]) > UMBRAL_COLIN)
    pares <- rbind(pares, data.frame(A=bandas_opt[i], B=bandas_opt[j],
                                     rho=round(mat_opt[i,j],3)))
}
eliminar <- c()
if (nrow(pares) > 0) {
  print(pares, row.names=FALSE)
  for (k in 1:nrow(pares)) {
    a <- pares$A[k]; b <- pares$B[k]
    if (a %in% eliminar || b %in% eliminar) next
    if (univ_opt[a] >= univ_opt[b]) eliminar <- c(eliminar, b) else eliminar <- c(eliminar, a)
  }
}
bandas_opt_dep <- setdiff(bandas_opt, eliminar)
cat(sprintf("\nÓptico depurado (%d de %d): %s\n",
            length(bandas_opt_dep), length(bandas_opt),
            paste(bandas_opt_dep, collapse=", ")))

# -----------------------------------------------------------------------------
# 2-3. Modelos: óptico y fusión
# -----------------------------------------------------------------------------
bandas_fusion <- c(bandas_sar_dep, bandas_opt_dep)

m_opt    <- correr_lofocv(dataset, bandas_opt_dep, "ÓPTICO (benchmark)")
m_fusion <- correr_lofocv(dataset, bandas_fusion,  "FUSIÓN SAR + Óptico")

# Recuperar SAR ya calculado
sar <- readRDS(RUTAS$rds_sar)
m_sar_global <- sar$metricas_global

# -----------------------------------------------------------------------------
# 4. Tabla comparativa final
# -----------------------------------------------------------------------------
cat("\n\n=============================================================\n")
cat("COMPARATIVA FINAL · SEVERIDAD (CBI_total)\n")
cat("=============================================================\n")
comp <- rbind(
  data.frame(Modelo="SAR (8 bandas dep.)",   m_sar_global[,c("r2","r_pearson","rmse","mae")]),
  data.frame(Modelo="Óptico (benchmark)",     m_opt$global[,c("r2","r_pearson","rmse","mae")]),
  data.frame(Modelo="Fusión SAR+Óptico",      m_fusion$global[,c("r2","r_pearson","rmse","mae")])
)
comp[,2:5] <- round(comp[,2:5], 3)
print(comp, row.names = FALSE)

cat("\n--- Lectura ---\n")
cat("Óptico vs SAR: cuánto supera el óptico al radar en graduar severidad.\n")
cat("Fusión vs Óptico: si la fusión supera al óptico solo, el SAR aporta\n")
cat("                  información estructural complementaria (Belenguer-Plomer 2021).\n")
cat("                  Si no lo supera, el óptico domina y el SAR no añade.\n")

# -----------------------------------------------------------------------------
# 5. Importancia de variables en el modelo de FUSIÓN
#    ¿El RF usa las bandas SAR o las ignora frente al óptico?
# -----------------------------------------------------------------------------
cat("\n\n=============================================================\n")
cat("IMPORTANCIA DE VARIABLES EN LA FUSIÓN\n")
cat("=============================================================\n")
set.seed(SEED)
form_fus <- as.formula(paste("CBI_total ~", paste(bandas_fusion, collapse=" + ")))
rf_fus <- ranger(form_fus, data = dataset, num.trees = N_ARBOLES,
                 importance = "permutation", seed = SEED)
imp_fus <- sort(rf_fus$variable.importance, decreasing = TRUE)

# Marcar cada variable como SAR u Óptico
tipo <- ifelse(names(imp_fus) %in% bandas_sar_dep, "SAR", "Óptico")
tabla_imp <- data.frame(
  Variable = names(imp_fus),
  Tipo = tipo,
  Importancia = round(as.numeric(imp_fus), 4)
)
print(tabla_imp, row.names = FALSE)

# Resumen agregado por tipo
cat("\n--- Importancia agregada por sensor ---\n")
imp_sar_total <- sum(imp_fus[names(imp_fus) %in% bandas_sar_dep])
imp_opt_total <- sum(imp_fus[names(imp_fus) %in% bandas_opt_dep])
imp_total <- imp_sar_total + imp_opt_total
cat(sprintf("  SAR    : %.4f (%.1f%% del total)\n",
            imp_sar_total, 100*imp_sar_total/imp_total))
cat(sprintf("  Óptico : %.4f (%.1f%% del total)\n",
            imp_opt_total, 100*imp_opt_total/imp_total))

# Posición de la primera variable SAR en el ranking
primera_sar <- which(tipo == "SAR")[1]
cat(sprintf("\n  La primera variable SAR aparece en la posición %d de %d.\n",
            primera_sar, length(imp_fus)))
cat("  Si las SAR están al fondo del ranking y suman poco %, se confirma\n")
cat("  cuantitativamente que el SAR no aporta sobre el óptico.\n")
saveRDS(list(comparativa=comp, optico=m_opt, fusion=m_fusion,
             bandas_opt_dep=bandas_opt_dep, bandas_fusion=bandas_fusion),
        file.path(RUTAS$dir_salida, "optico_fusion_severidad_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "optico_fusion_severidad_v3.rds")))

