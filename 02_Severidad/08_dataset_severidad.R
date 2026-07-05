# =============================================================================
# TFM SENTINEL-1 · DATASET DE SEVERIDAD (OBJETIVO 2)
# =============================================================================
# Construye el dataset para modelar CBI_total:
#   1. Carga 185 parcelas CBI de campo (8 incendios; descartados Cuevas,
#      SanBartolome, Ciperez).
#   2. Selecciona 50 puntos CBI=0 del muestreo fotointerpretado (quemado=0,
#      zona=fuera, a >=50 m del perímetro), repartidos proporcionalmente entre
#      incendios. Anclan el extremo bajo del gradiente de severidad.
#   3. Extrae, con buffer de 10 m (coherente con la celda de 20x20 m del CBI),
#      las 14 bandas SAR y las 11 ópticas en todos los registros.
#   4. Guarda el dataset listo para univariante + colinealidad + modelos.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
})

RUTAS <- list(
  csv_cbi     = "G:/TFM/_DATOS/CBI_incendios.csv",
  shp_puntos  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/puntos_fotointerpretacion_leon.shp",
  shp_perim   = "G:/TFM/_CARTOGRAFIA/PERIMETROS_INCENDIOS/Perimetros_CyL.shp",
  dir_s1      = "G:/TFM/_DATOS/EXPORT_S1_TFM",
  dir_s2      = "G:/TFM/_DATOS/EXPORT_S2_TFM",
  dir_salida  = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS"
)

INCENDIOS_EXCLUIR <- c("Cuevas", "SanBartolome", "Ciperez")
N_PUNTOS_CERO  <- 50
DIST_MIN_PERIM <- 50      # m mínimos al perímetro para puntos CBI=0
BUFFER_M       <- 10      # buffer de extracción (celda CBI ~20x20 m)
SEED           <- 42
set.seed(SEED)

bandas_sar <- c("Delta_VH_dB","Delta_VV_dB","Delta_SPAN_dB",
                "RBR_VH_Lin","RBR_VV_Lin","Delta_RFDI","Delta_RVI",
                "STD30_DeltaVH","STD70_DeltaVH","GLCM5_contrast","GLCM5_entropy",
                "GLCM5_variance","pre_VH_dB","pre_VV_dB")
bandas_opt <- c("dNBR","RBR","RdNBR","dNDMI","dBAIS2","dNBRplus","dNDVI",
                "pre_NBR","pre_NDMI","pre_BAIS2","pre_NDVI")

cat("=============================================================\n")
cat("DATASET DE SEVERIDAD (OBJETIVO 2)\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# 1. Parcelas CBI
# -----------------------------------------------------------------------------
cbi <- read.csv(RUTAS$csv_cbi)
cbi <- cbi[!cbi$Incendio %in% INCENDIOS_EXCLUIR, ]
cat(sprintf("Parcelas CBI (8 incendios): %d\n", nrow(cbi)))

parcelas <- data.frame(
  origen   = "campo",
  Incendio = cbi$Incendio,
  x_utm    = cbi$X,
  y_utm    = cbi$Y,
  CBI_total = cbi$CBI_total
)

# -----------------------------------------------------------------------------
# 2. Puntos CBI=0 (fuera del perímetro, alejados, repartidos por incendio)
# -----------------------------------------------------------------------------
pts <- st_read(RUTAS$shp_puntos, quiet = TRUE)
if (st_crs(pts)$epsg != 25830) pts <- st_transform(pts, 25830)
perim <- st_read(RUTAS$shp_perim, quiet = TRUE)
if (st_crs(perim)$epsg != 25830) perim <- st_transform(perim, 25830)

# Candidatos: quemado=0, zona=fuera
cand <- pts[!is.na(pts$quemado) & pts$quemado == 0 & pts$zona == "fuera", ]
cand <- cand[!cand$Incendio %in% INCENDIOS_EXCLUIR, ]

# Distancia al perímetro de su incendio >= DIST_MIN_PERIM
cand$dist_perim <- NA_real_
for (inc in unique(cand$Incendio)) {
  per_inc <- st_union(perim[perim$COMMUNE == inc, ])
  idx <- which(cand$Incendio == inc)
  if (length(idx) > 0 && length(per_inc) > 0) {
    cand$dist_perim[idx] <- as.numeric(st_distance(cand[idx, ], per_inc))
  }
}
cand <- cand[!is.na(cand$dist_perim) & cand$dist_perim >= DIST_MIN_PERIM, ]
cat(sprintf("Candidatos CBI=0 (fuera, >=%dm del perímetro): %d\n",
            DIST_MIN_PERIM, nrow(cand)))

# Reparto proporcional al número de candidatos por incendio
tab <- table(cand$Incendio)
prop <- as.numeric(tab) / sum(tab)
n_por_inc <- round(prop * N_PUNTOS_CERO)
names(n_por_inc) <- names(tab)
# Ajuste para cuadrar exactamente N_PUNTOS_CERO
while (sum(n_por_inc) != N_PUNTOS_CERO) {
  if (sum(n_por_inc) < N_PUNTOS_CERO) {
    i <- which.max(prop); n_por_inc[i] <- n_por_inc[i] + 1
  } else {
    i <- which.max(n_por_inc); n_por_inc[i] <- n_por_inc[i] - 1
  }
}

cat("\n--- Reparto de puntos CBI=0 por incendio ---\n")
sel_idx <- c()
for (inc in names(n_por_inc)) {
  idx_inc <- which(cand$Incendio == inc)
  n_pedir <- min(n_por_inc[inc], length(idx_inc))
  if (n_pedir < n_por_inc[inc]) {
    cat(sprintf("  [aviso] %s: solo %d candidatos (<%d pedidos)\n",
                inc, length(idx_inc), n_por_inc[inc]))
  }
  sel_idx <- c(sel_idx, sample(idx_inc, n_pedir))
  cat(sprintf("  %-12s: %d puntos\n", inc, n_pedir))
}
cand_sel <- cand[sel_idx, ]
cc <- st_coordinates(cand_sel)

puntos_cero <- data.frame(
  origen   = "cero",
  Incendio = cand_sel$Incendio,
  x_utm    = cc[,1],
  y_utm    = cc[,2],
  CBI_total = 0
)
cat(sprintf("\nPuntos CBI=0 seleccionados: %d\n", nrow(puntos_cero)))

# -----------------------------------------------------------------------------
# 3. Unir parcelas + puntos CBI=0
# -----------------------------------------------------------------------------
dataset <- bind_rows(parcelas, puntos_cero)
cat(sprintf("Dataset total: %d registros (%d campo + %d cero)\n",
            nrow(dataset), sum(dataset$origen=="campo"),
            sum(dataset$origen=="cero")))

# Convertir a sf
ds_sf <- st_as_sf(dataset, coords = c("x_utm","y_utm"), crs = 25830, remove = FALSE)

# -----------------------------------------------------------------------------
# 4. Extracción por buffer 10 m de SAR y ópticas
# -----------------------------------------------------------------------------
extraer_bandas <- function(dir_tif, bandas, etiqueta) {
  tifs <- list.files(dir_tif, pattern = "\\.tif$", full.names = TRUE)
  res_list <- list()
  for (inc in unique(ds_sf$Incendio)) {
    tp <- tifs[grepl(paste0(inc, "\\.tif"), basename(tifs))]
    if (length(tp) == 0) {
      cat(sprintf("  [%s] sin TIF %s; se salta.\n", inc, etiqueta)); next
    }
    r <- rast(tp[1])
    bandas_ok <- intersect(bandas, names(r))
    r <- r[[bandas_ok]]
    pts_inc <- ds_sf[ds_sf$Incendio == inc, ]
    buf <- st_buffer(pts_inc, BUFFER_M)
    v <- terra::extract(r, vect(buf), fun = mean, na.rm = TRUE, ID = FALSE)
    v$row_id <- which(ds_sf$Incendio == inc)
    res_list[[inc]] <- v
  }
  do.call(rbind, res_list)
}

cat("\nExtrayendo bandas SAR (buffer 10 m)...\n")
ext_sar <- extraer_bandas(RUTAS$dir_s1, bandas_sar, "SAR")
cat("Extrayendo bandas ópticas (buffer 10 m)...\n")
ext_opt <- extraer_bandas(RUTAS$dir_s2, bandas_opt, "OPT")

# Reensamblar en el orden original
ext_sar <- ext_sar[order(ext_sar$row_id), ]
ext_opt <- ext_opt[order(ext_opt$row_id), ]
dataset_final <- cbind(
  dataset,
  ext_sar[, bandas_sar, drop = FALSE],
  ext_opt[, bandas_opt, drop = FALSE]
)

# Eliminar registros con NA en predictores
ok <- complete.cases(dataset_final[, c(bandas_sar, bandas_opt)])
n_descart <- sum(!ok)
dataset_final <- dataset_final[ok, ]
cat(sprintf("\nRegistros descartados por NA en predictores: %d\n", n_descart))
cat(sprintf("Dataset final: %d registros\n", nrow(dataset_final)))

# -----------------------------------------------------------------------------
# Resumen y guardado
# -----------------------------------------------------------------------------
cat("\n--- Distribución por incendio ---\n")
res <- dataset_final %>%
  group_by(Incendio) %>%
  summarise(n_total = n(),
            n_campo = sum(origen=="campo"),
            n_cero  = sum(origen=="cero"),
            cbi_medio = round(mean(CBI_total),2),
            .groups="drop")
print(as.data.frame(res), row.names = FALSE)

saveRDS(list(
  dataset = dataset_final,
  bandas_sar = bandas_sar,
  bandas_opt = bandas_opt,
  buffer_m = BUFFER_M,
  n_puntos_cero = N_PUNTOS_CERO
), file.path(RUTAS$dir_salida, "dataset_severidad_v3.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "dataset_severidad_v3.rds")))
cat("Dataset de severidad completado.\n")

