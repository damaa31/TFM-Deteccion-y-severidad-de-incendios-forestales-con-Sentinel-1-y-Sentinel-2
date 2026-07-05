# =============================================================================
# TFM SENTINEL-1 · MÓDULO 2 (v3): EXTRACCIÓN SAR EN PUNTOS FOTOINTERPRETADOS
# =============================================================================
# Construye el dataset de modelado de detección:
#   - Carga el shapefile con los puntos fotointerpretados (campo 'quemado').
#   - Para cada punto, extrae el valor de las 14 bandas Sentinel-1 en su píxel.
#   - Genera la tabla final lista para modelar (puntos x predictores + clase).
#
# Notas:
#   - Extracción puntual (un único píxel por punto). Coherente con la unidad
#     de fotointerpretación.
#   - Se eliminan los puntos dudosos (NA en 'quemado').
#   - Se eliminan los puntos sin valor SAR válido (raros, en bordes del TIF).
#
# Salida: dataset_modelado_v3.rds
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
})

RUTAS <- list(
  shp_puntos = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/puntos_fotointerpretacion_leon.shp",
  dir_s1     = "G:/TFM/_DATOS/EXPORT_S1_TFM/",
  dir_salida = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/"
)

# Nombre del campo con la etiqueta fotointerpretada (0/1, NA dudoso)
CAMPO_CLASE <- "quemado"

# Bandas SAR a extraer (las 14 estándar del estudio)
bandas_sar <- c("Delta_VH_dB","Delta_VV_dB","Delta_SPAN_dB",
                "RBR_VH_Lin","RBR_VV_Lin","Delta_RFDI","Delta_RVI",
                "STD30_DeltaVH","STD70_DeltaVH","GLCM5_contrast","GLCM5_entropy",
                "GLCM5_variance","pre_VH_dB","pre_VV_dB")

cat("=============================================================\n")
cat("MÓDULO 2 (v3): EXTRACCIÓN SAR EN PUNTOS FOTOINTERPRETADOS\n")
cat("=============================================================\n")

# -----------------------------------------------------------------------------
# Cargar puntos
# -----------------------------------------------------------------------------
pts <- st_read(RUTAS$shp_puntos, quiet = TRUE)
if (st_crs(pts)$epsg != 25830) pts <- st_transform(pts, 25830)
cat(sprintf("Puntos cargados: %d\n", nrow(pts)))

if (!CAMPO_CLASE %in% names(pts)) {
  stop(sprintf("No se encuentra el campo '%s' en el shapefile.", CAMPO_CLASE))
}

# Filtrar dudosos (NA en la etiqueta)
n_na <- sum(is.na(pts[[CAMPO_CLASE]]))
pts <- pts[!is.na(pts[[CAMPO_CLASE]]), ]
cat(sprintf("Puntos dudosos descartados: %d | Puntos válidos: %d\n",
            n_na, nrow(pts)))

# -----------------------------------------------------------------------------
# TIFs disponibles
# -----------------------------------------------------------------------------
tifs_s1 <- list.files(RUTAS$dir_s1, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("TIFs Sentinel-1 disponibles: %d\n\n", length(tifs_s1)))

# -----------------------------------------------------------------------------
# Extracción por incendio
# -----------------------------------------------------------------------------
lista_extracciones <- list()

for (id in unique(pts$Incendio)) {
  cat(sprintf("--- [%s] ---\n", id))
  
  # Localizar TIF del incendio (tolerante a .tif y .tif.tif)
  patron <- paste0(id, "\\.tif")
  tif_path <- tifs_s1[grepl(patron, basename(tifs_s1))]
  if (length(tif_path) == 0) {
    cat("  Sin TIFF SAR; se salta.\n"); next
  }
  
  # Cargar el TIF y quedarnos con las 14 bandas (en el orden esperado)
  r <- rast(tif_path[1])
  bandas_disponibles <- intersect(bandas_sar, names(r))
  faltan <- setdiff(bandas_sar, names(r))
  if (length(faltan) > 0) {
    cat(sprintf("  [aviso] faltan bandas: %s\n", paste(faltan, collapse = ", ")))
  }
  r <- r[[bandas_disponibles]]
  
  # Puntos del incendio
  pts_inc <- pts[pts$Incendio == id, ]
  cat(sprintf("  Puntos a extraer: %d\n", nrow(pts_inc)))
  
  # Extracción puntual (un único píxel por punto)
  vals <- terra::extract(r, vect(pts_inc), ID = FALSE)
  
  # Ensamblar
  df_inc <- data.frame(
    Incendio   = id,
    punto_id   = pts_inc$punto_id,
    x_utm      = pts_inc$x_utm,
    y_utm      = pts_inc$y_utm,
    quemado    = pts_inc[[CAMPO_CLASE]]
  )
  df_inc <- cbind(df_inc, vals)
  lista_extracciones[[id]] <- df_inc
  
  cat(sprintf("  Filas extraídas: %d\n", nrow(df_inc)))
}

# -----------------------------------------------------------------------------
# Consolidar y filtrar puntos sin SAR válido
# -----------------------------------------------------------------------------
dataset <- bind_rows(lista_extracciones)
cat(sprintf("\nDataset consolidado: %d filas\n", nrow(dataset)))

# Eliminar puntos con NA en alguna banda SAR
ok <- complete.cases(dataset[, bandas_sar])
n_descart <- sum(!ok)
dataset <- dataset[ok, ]
cat(sprintf("Puntos descartados por SAR no válido: %d\n", n_descart))
cat(sprintf("Dataset final: %d filas\n", nrow(dataset)))

# -----------------------------------------------------------------------------
# Resumen
# -----------------------------------------------------------------------------
cat("\n--- Distribución de clases por incendio (dataset final) ---\n")
res <- dataset %>%
  group_by(Incendio) %>%
  summarise(
    n_total   = n(),
    n_quemado = sum(quemado == 1),
    n_no_quem = sum(quemado == 0),
    pct_quem  = round(100 * n_quemado / n_total, 1),
    .groups   = "drop"
  )
print(as.data.frame(res), row.names = FALSE)

cat(sprintf("\nTotal: %d puntos (%d quemado, %d no quemado)\n",
            nrow(dataset),
            sum(dataset$quemado == 1),
            sum(dataset$quemado == 0)))

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
saveRDS(list(
  dataset = dataset,
  bandas_sar = bandas_sar,
  resumen_por_incendio = res
), file.path(RUTAS$dir_salida, "dataset_modelado_v3.rds"))

cat(sprintf("\nGuardado: %s\n",
            file.path(RUTAS$dir_salida, "dataset_modelado_v3.rds")))
cat("Módulo 2 v3 completado.\n")