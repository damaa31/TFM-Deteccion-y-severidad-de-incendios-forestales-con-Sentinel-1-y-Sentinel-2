# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · OPTICO POR ECOSISTEMA (arbolado vs matorral)
# =============================================================================
# Replica EXACTAMENTE la metodologia aplicada al SAR (script 12 y 13) pero con
# el optico, para comparar sensores por cubierta en igualdad de condiciones:
#   - Optico depurado por colinealidad (Spearman |rho|>0.90), criterio identico.
#   - Univariante Spearman por ecosistema (solo parcelas de campo).
#   - Modelos RF k-fold por ecosistema (campo + puntos CBI=0 clasificados).
#   - Comprobacion LOFOCV de robustez.
#
# Pregunta: el optico, ?es bueno en ambos ecosistemas, o tambien depende de la
# cubierta? Si el optico es bueno en ambos y el SAR solo en matorral, el
# matorral es el terreno donde el SAR se aproxima al optico.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(ranger); library(dplyr)
})

RUTAS <- list(
  csv_cbi     = "G:/TFM/_DATOS/CBI_incendios.csv",
  rds_dataset = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/dataset_severidad_v3.rds",
  mf_raster   = "G:/TFM/_CARTOGRAFIA/FORESTAL/MF_raster.tif",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)
INCENDIOS_EXCLUIR <- c("Cuevas","SanBartolome","Ciperez")
N_ARBOLES <- 500
K_FOLDS   <- 5
SEED <- 42
UMBRAL_COLIN <- 0.90
MF_COD <- c("1"="arbolado", "2"="matorral")

rds  <- readRDS(RUTAS$rds_dataset)
dataset    <- rds$dataset
bandas_opt <- rds$bandas_opt

cat("=============================================================\n")
cat("OPTICO POR ECOSISTEMA (arbolado vs matorral)\n")
cat("=============================================================\n")

# --- Asignar ecosistema (identico a script 12) ---
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

# --- Depuracion colinealidad del optico (criterio identico al SAR) ---
cat("\n--- Depuracion colinealidad OPTICO (Spearman |rho|>0.90) ---\n")
univ_opt_full <- sapply(bandas_opt, function(b)
  abs(cor(dataset[[b]], dataset$CBI_total, method="spearman")))
mat_opt <- cor(dataset[, bandas_opt], method="spearman", use="complete.obs")
eliminar <- c()
for (i in 1:(length(bandas_opt)-1)) for (j in (i+1):length(bandas_opt)) {
  if (abs(mat_opt[i,j]) > UMBRAL_COLIN) {
    a <- bandas_opt[i]; b <- bandas_opt[j]
    if (a %in% eliminar || b %in% eliminar) next
    if (univ_opt_full[a] >= univ_opt_full[b]) eliminar <- c(eliminar,b) else eliminar <- c(eliminar,a)
  }
}
bandas_opt_dep <- setdiff(bandas_opt, eliminar)
cat(sprintf("Optico depurado (%d de %d): %s\n",
            length(bandas_opt_dep), length(bandas_opt),
            paste(bandas_opt_dep, collapse=", ")))

calc_metricas <- function(obs, pred) {
  ss_res <- sum((obs-pred)^2); ss_tot <- sum((obs-mean(obs))^2)
  data.frame(r2=1-ss_res/ss_tot, r_pearson=cor(obs,pred),
             rmse=sqrt(mean((obs-pred)^2)), mae=mean(abs(obs-pred)))
}

# =============================================================================
# C - UNIVARIANTE OPTICO POR ECOSISTEMA (solo campo)
# =============================================================================
cat("\n=============================================================\n")
cat("UNIVARIANTE OPTICO POR ECOSISTEMA (Spearman, solo campo)\n")
cat("=============================================================\n")
dat_campo <- dataset[campo & !is.na(dataset$eco), ]
univ_eco <- data.frame()
for (b in bandas_opt) {
  rho_a <- cor(dat_campo[[b]][dat_campo$eco=="arbolado"],
               dat_campo$CBI_total[dat_campo$eco=="arbolado"], method="spearman")
  rho_m <- cor(dat_campo[[b]][dat_campo$eco=="matorral"],
               dat_campo$CBI_total[dat_campo$eco=="matorral"], method="spearman")
  univ_eco <- rbind(univ_eco, data.frame(
    Banda=b, rho_arbolado=round(rho_a,3), rho_matorral=round(rho_m,3),
    dif_abs=round(abs(rho_m)-abs(rho_a),3)))
}
univ_eco <- univ_eco[order(-abs(univ_eco$rho_matorral)), ]
print(univ_eco, row.names=FALSE)
cat(sprintf("\n  |rho| medio arbolado: %.3f | matorral: %.3f\n",
            mean(abs(univ_eco$rho_arbolado)), mean(abs(univ_eco$rho_matorral))))

# =============================================================================
# B - MODELOS OPTICO k-fold + LOFOCV por ecosistema
# =============================================================================
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
  m <- calc_metricas(pred_all$obs, pred_all$pred); m$grupo<-etiqueta; m$n<-nrow(dat); m
}
lofocv <- function(dat, bandas, etiqueta) {
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  pred_all <- data.frame()
  for (inc in unique(dat$Incendio)) {
    train <- dat[dat$Incendio!=inc,]; test <- dat[dat$Incendio==inc,]
    if (nrow(test)==0) next
    set.seed(SEED)
    rf <- ranger(form, data=train, num.trees=N_ARBOLES, seed=SEED)
    pred <- pmin(pmax(predict(rf,test)$predictions,0),3)
    pred_all <- rbind(pred_all, data.frame(obs=test$CBI_total, pred=pred))
  }
  m <- calc_metricas(pred_all$obs, pred_all$pred); m$grupo<-etiqueta; m$n<-nrow(pred_all); m
}

dat_eco <- dataset[!is.na(dataset$eco), ]
cat("\n=============================================================\n")
cat("MODELOS OPTICO POR ECOSISTEMA\n")
cat("=============================================================\n")
kf <- rbind(
  modelo_kfold(dat_eco[dat_eco$eco=="arbolado",], bandas_opt_dep, "Optico arbolado"),
  modelo_kfold(dat_eco[dat_eco$eco=="matorral",], bandas_opt_dep, "Optico matorral")
)
lo <- rbind(
  lofocv(dat_eco[dat_eco$eco=="arbolado",], bandas_opt_dep, "Optico arbolado"),
  lofocv(dat_eco[dat_eco$eco=="matorral",], bandas_opt_dep, "Optico matorral")
)
cat("\n--- k-fold ---\n")
print(kf[,c("grupo","n","r2","r_pearson","rmse","mae")] %>%
        mutate(across(where(is.numeric),~round(.,3))), row.names=FALSE)
cat("\n--- LOFOCV (robustez) ---\n")
print(lo[,c("grupo","n","r2","r_pearson","rmse","mae")] %>%
        mutate(across(where(is.numeric),~round(.,3))), row.names=FALSE)

cat("\n=============================================================\n")
cat("COMPARATIVA SENSORES POR ECOSISTEMA (LOFOCV)\n")
cat("=============================================================\n")
cat("  SAR    : arbolado R2=-0.025 | matorral R2=0.426\n")
cat(sprintf("  Optico : arbolado R2=%.3f | matorral R2=%.3f\n",
            lo$r2[lo$grupo=="Optico arbolado"], lo$r2[lo$grupo=="Optico matorral"]))

saveRDS(list(univariante_opt_eco=univ_eco, kfold=kf, lofocv=lo,
             bandas_opt_dep=bandas_opt_dep),
        file.path(RUTAS$dir_salida, "optico_por_ecosistema_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "optico_por_ecosistema_v3.rds")))