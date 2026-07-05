# =============================================================================
# TFM SENTINEL-1 · PUNTOS ALEATORIOS PARA FOTOINTERPRETACIÓN
# =============================================================================
# Diseño:
#   - 150 puntos por incendio (Padilla et al. 2015; Olofsson et al. 2014).
#   - Aleatorios sobre el bbox del TIF, clasificados a posteriori en
#     dentro/fuera del perímetro vía st_intersects (estable).
#   - Mínimo 40 dentro / 40 fuera; el resto proporcional al reparto natural.
#   - Distancia mínima entre puntos: 50 m.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
})

RUTAS <- list(
  shp_perimetros = "G:/TFM/_CARTOGRAFIA/PERIMETROS_INCENDIOS/Perimetros_CyL.shp",
  dir_s1         = "G:/TFM/_DATOS/EXPORT_S1_TFM/",
  dir_salida     = "G:/TFM/_DATOS/_ANALISIS/_RESULTADOS/"
)
CAMPO_ID <- "COMMUNE"
INCENDIOS_EXCLUIR <- c("Cuevas", "SanBartolome")

N_POR_INCENDIO <- 150
MIN_DENTRO     <- 40
MIN_FUERA      <- 40
N_MIN_INCENDIO <- 100   # mínimo absoluto por incendio
DIST_MIN_M     <- 50
BASE_N         <- 150     # base de referencia por incendio (TIF mediano)
AREA_BASE_HA   <- NULL    # se calculará al cargar los TIFs
SEED           <- 42
set.seed(SEED)

dir.create(RUTAS$dir_salida, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================\n")
cat("PUNTOS ALEATORIOS PARA FOTOINTERPRETACIÓN\n")
cat("=============================================================\n")

perim <- st_read(RUTAS$shp_perimetros, quiet = TRUE)
if (st_crs(perim)$epsg != 25830) perim <- st_transform(perim, 25830)
perim <- perim[!(perim[[CAMPO_ID]] %in% INCENDIOS_EXCLUIR), ]
cat(sprintf("Incendios: %s\n\n", paste(perim[[CAMPO_ID]], collapse = ", ")))

tifs_s1 <- list.files(RUTAS$dir_s1, pattern = "\\.tif$", full.names = TRUE)

# -----------------------------------------------------------------------------
# Muestreo aleatorio uniforme dentro de un bbox, con distancia mínima.
# Trabaja con coordenadas crudas, sin sf, para máxima estabilidad.
# -----------------------------------------------------------------------------
muestrear_bbox <- function(xmin, xmax, ymin, ymax, n, dist_min, max_iter = 50) {
  acum <- matrix(nrow = 0, ncol = 2)
  iter <- 0
  while (nrow(acum) < n && iter < max_iter) {
    n_pedir <- (n - nrow(acum)) * 5
    x <- runif(n_pedir, xmin, xmax)
    y <- runif(n_pedir, ymin, ymax)
    for (i in seq_len(n_pedir)) {
      p <- c(x[i], y[i])
      if (nrow(acum) == 0) {
        acum <- rbind(acum, p)
      } else {
        d <- sqrt((acum[,1] - p[1])^2 + (acum[,2] - p[2])^2)
        if (all(d >= dist_min)) acum <- rbind(acum, p)
      }
      if (nrow(acum) >= n) break
    }
    iter <- iter + 1
  }
  acum
}

# -----------------------------------------------------------------------------
# Procesar cada incendio
# -----------------------------------------------------------------------------
lista_pts <- list()
resumen <- data.frame()
# Calcular el área del TIF mediano para escalar proporcionalmente
areas_tifs <- sapply(seq_len(nrow(perim)), function(k) {
  id <- as.character(perim[[CAMPO_ID]][k])
  patron <- paste0(id, "\\.tif")
  tp <- tifs_s1[grepl(patron, basename(tifs_s1))]
  if (length(tp) == 0) return(NA_real_)
  r <- rast(tp[1])[[1]]
  bb <- ext(r)
  (bb$xmax - bb$xmin) * (bb$ymax - bb$ymin) / 10000  # ha
})
AREA_BASE_HA <- median(areas_tifs, na.rm = TRUE)
cat(sprintf("Área TIF de referencia (mediana): %.0f ha\n\n", AREA_BASE_HA))

for (k in seq_len(nrow(perim))) {
  id <- as.character(perim[[CAMPO_ID]][k])
  cat(sprintf("--- [%s] ---\n", id))
  
  patron <- paste0(id, "\\.tif")
  tif_path <- tifs_s1[grepl(patron, basename(tifs_s1))]
  if (length(tif_path) == 0) {
    cat("  Sin TIFF SAR; se salta.\n"); next
  }
  r <- rast(tif_path[1])[[1]]
  bb <- ext(r)
  per_inc <- st_geometry(perim[k, ])
  
  # Sobre-muestrear (5x el N total) y luego clasificar y reequilibrar.
  area_tif_ha_temp <- (bb$xmax - bb$xmin) * (bb$ymax - bb$ymin) / 10000
  N_SOBRE <- max(750, round(BASE_N * area_tif_ha_temp / AREA_BASE_HA) * 5)
  coords <- muestrear_bbox(bb$xmin, bb$xmax, bb$ymin, bb$ymax,
                           N_SOBRE, DIST_MIN_M)
  cat(sprintf("  Sobre-muestreo: %d puntos\n", nrow(coords)))
  
  # Convertir a sf y clasificar dentro/fuera vía st_intersects (estable)
  pts_all <- st_as_sf(data.frame(x = coords[,1], y = coords[,2]),
                      coords = c("x","y"), crs = 25830, remove = FALSE)
  hits <- lengths(st_intersects(pts_all, per_inc)) > 0
  pts_all$zona <- ifelse(hits, "dentro", "fuera")
  
  n_dentro_disp <- sum(pts_all$zona == "dentro")
  n_fuera_disp  <- sum(pts_all$zona == "fuera")
  cat(sprintf("  Disponibles: dentro=%d | fuera=%d\n",
              n_dentro_disp, n_fuera_disp))
  
  # N total proporcional al área del TIF respecto al TIF mediano,
  # con suelo mínimo para garantizar potencia estadística por incendio
  area_tif_ha <- (bb$xmax - bb$xmin) * (bb$ymax - bb$ymin) / 10000
  n_total_prop <- round(BASE_N * area_tif_ha / AREA_BASE_HA)
  n_total <- max(N_MIN_INCENDIO, n_total_prop)
  
  # Reparto 50/50 dentro/fuera del perímetro
  n_dentro <- min(round(n_total / 2), n_dentro_disp)
  n_fuera  <- min(n_total - n_dentro, n_fuera_disp)
  cat(sprintf("  N proporcional: total=%d (área TIF %.0f ha vs base %.0f ha) | dentro=%d | fuera=%d\n",
              n_total, area_tif_ha, AREA_BASE_HA, n_dentro, n_fuera))
  
  # Seleccionar aleatoriamente N de cada estrato
  idx_d <- which(pts_all$zona == "dentro")
  idx_f <- which(pts_all$zona == "fuera")
  sel_d <- sample(idx_d, n_dentro)
  sel_f <- sample(idx_f, n_fuera)
  pts_sf <- pts_all[c(sel_d, sel_f), ]
  
  pts_sf$Incendio  <- id
  pts_sf$punto_id  <- sprintf("%s_%04d", id, seq_len(nrow(pts_sf)))
  pts_sf$y_quemado <- NA_integer_
  cc <- st_coordinates(pts_sf)
  pts_sf$x_utm <- cc[,1]; pts_sf$y_utm <- cc[,2]
  pts_sf$x <- NULL; pts_sf$y <- NULL   # quitar las originales redundantes
  
  lista_pts[[id]] <- pts_sf
  resumen <- rbind(resumen, data.frame(
    Incendio = id, n_total = nrow(pts_sf),
    n_dentro = sum(pts_sf$zona == "dentro"),
    n_fuera  = sum(pts_sf$zona == "fuera")))
}

# -----------------------------------------------------------------------------
# Guardado
# -----------------------------------------------------------------------------
if (length(lista_pts) == 0) stop("No se generaron puntos.")

todos <- do.call(rbind, lista_pts)
cat("\n=== RESUMEN ===\n"); print(resumen, row.names = FALSE)
cat(sprintf("\nTotal puntos: %d\n", nrow(todos)))

out_shp <- file.path(RUTAS$dir_salida, "puntos_fotointerpretacion.shp")
st_write(todos, out_shp, delete_layer = TRUE, quiet = TRUE)
cat(sprintf("Shapefile: %s\n", out_shp))

out_csv <- file.path(RUTAS$dir_salida, "puntos_fotointerpretacion.csv")
write.csv(st_drop_geometry(todos), out_csv, row.names = FALSE)
cat(sprintf("CSV: %s\n\n", out_csv))

cat("Abre el shapefile en QGIS sobre PNOA o S2 RGB post-fuego.\n")
cat("Rellena la columna y_quemado: 1 = quemado, 0 = no quemado, NA = dudoso.\n")