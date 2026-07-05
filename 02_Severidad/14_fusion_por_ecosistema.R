# =============================================================================
# TFM SENTINEL-1 · OBJETIVO 2 · FUSION SAR+OPTICO POR ECOSISTEMA
# =============================================================================
# Cierra la cuestion de complementariedad: ?aporta el SAR a la fusion
# especificamente en MATORRAL, donde el SAR es competente (R2=0.426) y el optico
# algo mas debil (R2=0.666)? Replica la metodologia (depuracion identica,
# k-fold + LOFOCV por ecosistema) y compara optico solo vs fusion, por cubierta.
# Incluye importancia de variables por sensor en la fusion de matorral.
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

rds  <- readRDS(RUTAS$rds_dataset)
coln <- readRDS(RUTAS$rds_colin)
opt_eco <- readRDS(RUTAS$rds_opt_eco)
dataset    <- rds$dataset
bandas_sar_dep <- coln$bandas_depuradas
bandas_opt_dep <- opt_eco$bandas_opt_dep
bandas_fusion  <- c(bandas_sar_dep, bandas_opt_dep)

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
dat_eco <- dataset[!is.na(dataset$eco), ]

calc_metricas <- function(obs, pred) {
  ss_res <- sum((obs-pred)^2); ss_tot <- sum((obs-mean(obs))^2)
  data.frame(r2=1-ss_res/ss_tot, r_pearson=cor(obs,pred),
             rmse=sqrt(mean((obs-pred)^2)), mae=mean(abs(obs-pred)))
}
modelo_kfold <- function(dat, bandas) {
  set.seed(SEED)
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  folds <- sample(rep(1:K_FOLDS, length.out=nrow(dat)))
  pa <- data.frame()
  for (k in 1:K_FOLDS) {
    tr <- dat[folds!=k,]; te <- dat[folds==k,]
    rf <- ranger(form, data=tr, num.trees=N_ARBOLES, seed=SEED)
    pa <- rbind(pa, data.frame(obs=te$CBI_total,
                               pred=pmin(pmax(predict(rf,te)$predictions,0),3)))
  }
  calc_metricas(pa$obs, pa$pred)
}
lofocv <- function(dat, bandas) {
  form <- as.formula(paste("CBI_total ~", paste(bandas, collapse=" + ")))
  pa <- data.frame()
  for (inc in unique(dat$Incendio)) {
    tr <- dat[dat$Incendio!=inc,]; te <- dat[dat$Incendio==inc,]
    if (nrow(te)==0) next
    set.seed(SEED)
    rf <- ranger(form, data=tr, num.trees=N_ARBOLES, seed=SEED)
    pa <- rbind(pa, data.frame(obs=te$CBI_total,
                               pred=pmin(pmax(predict(rf,te)$predictions,0),3)))
  }
  calc_metricas(pa$obs, pa$pred)
}

cat("=============================================================\n")
cat("FUSION SAR+OPTICO POR ECOSISTEMA\n")
cat("=============================================================\n")

# --- Tabla comparativa: optico solo vs fusion, por ecosistema y validacion ---
resultados <- data.frame()
for (e in c("arbolado","matorral")) {
  d <- dat_eco[dat_eco$eco==e,]
  for (modelo in c("Optico","Fusion")) {
    bandas <- if (modelo=="Optico") bandas_opt_dep else bandas_fusion
    kf <- modelo_kfold(d, bandas); lo <- lofocv(d, bandas)
    resultados <- rbind(resultados, data.frame(
      Ecosistema=e, Modelo=modelo, n=nrow(d),
      R2_kfold=round(kf$r2,3), R2_lofocv=round(lo$r2,3),
      r_lofocv=round(lo$r_pearson,3), RMSE_lofocv=round(lo$rmse,3)))
  }
}
print(resultados, row.names=FALSE)

# --- Importancia por sensor en la FUSION de matorral ---
cat("\n--- Importancia por sensor en FUSION · MATORRAL ---\n")
set.seed(SEED)
d_mat <- dat_eco[dat_eco$eco=="matorral",]
form_f <- as.formula(paste("CBI_total ~", paste(bandas_fusion, collapse=" + ")))
rf_mat <- ranger(form_f, data=d_mat, num.trees=N_ARBOLES,
                 importance="permutation", seed=SEED)
imp <- sort(rf_mat$variable.importance, decreasing=TRUE)
tipo <- ifelse(names(imp) %in% bandas_sar_dep, "SAR", "Optico")
imp_sar <- sum(imp[names(imp) %in% bandas_sar_dep])
imp_opt <- sum(imp[names(imp) %in% bandas_opt_dep])
tot <- imp_sar + imp_opt
cat(sprintf("  SAR    : %.1f%% del total\n", 100*imp_sar/tot))
cat(sprintf("  Optico : %.1f%% del total\n", 100*imp_opt/tot))
cat(sprintf("  (En la fusion GLOBAL el SAR aportaba 5.3%%; aqui en matorral?)\n"))
primera_sar <- which(tipo=="SAR")[1]
cat(sprintf("  Primera variable SAR en posicion %d de %d.\n",
            primera_sar, length(imp)))

cat("\n--- Lectura ---\n")
cat("Si la fusion supera al optico SOLO en matorral, el SAR aporta justo donde\n")
cat("el optico es mas debil. Si no lo supera ni en matorral, el optico domina\n")
cat("en todo el rango y el SAR no anade ni en su mejor terreno.\n")

saveRDS(list(comparativa=resultados, importancia_matorral=imp),
        file.path(RUTAS$dir_salida, "fusion_por_ecosistema_v3.rds"))
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "fusion_por_ecosistema_v3.rds")))