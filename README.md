# TFM-Deteccion-y-severidad-de-incendios-forestales-con-Sentinel1-y-Sentinel2
Código del Trabajo Fin de Máster que evalúa el radar de banda C (Sentinel-1)  frente al sensor óptico (Sentinel-2) para la detección del área quemada y la  caracterización de la severidad del fuego

## Estructura del repositorio

Los scripts están numerados según el orden de ejecución del flujo de trabajo.

### 00_preprocesamiento_GEE — Google Earth Engine (JavaScript)
- `01_sentinel1_preproceso_stack.js` — Preprocesamiento y stack de Sentinel-1
- `02_sentinel2_stack_severidad.js` — Composición y stack de Sentinel-2

### 01_deteccion — Objetivo 1 (R)
- `03_puntos_aleatorios_fotointerpretacion.R` — Generación de puntos para fotointerpretación
- `04_extraccion_sar_puntos.R` — Extracción de variables SAR en los puntos
- `05_boruta_seleccion_variables.R` — Selección de variables con Boruta
- `06_modelo_deteccion_sar.R` — Modelo Random Forest de detección (SAR)
- `07_modelo_deteccion_optico_control.R` — Modelo óptico de control

### 02_severidad — Objetivo 2 (R)
- `08_dataset_severidad.R` — Construcción del dataset de severidad (CBI)
- `09_univariante_colinealidad_sar.R` — Análisis univariante y depuración de colinealidad
- `10_modelo_severidad_sar.R` — Modelo SAR de severidad continua
- `11_severidad_optico_fusion.R` — Modelos óptico y de fusión
- `12_sar_por_ecosistema.R` — Severidad SAR por ecosistema (arbolado/matorral)
- `13_optico_por_ecosistema.R` — Severidad óptica por ecosistema
- `14_fusion_por_ecosistema.R` — Fusión por ecosistema
- `15_clasificacion_severidad_categorica.R` — Clasificación categórica del CBI

### 03_cartografia — Mapas finales (R)
- `16_cartografia_deteccion.R` — Cartografía de detección (Objetivo 1)
- `17_cartografia_severidad.R` — Cartografía de severidad (Objetivo 2)
- `18_cartografia_severidad_fusion.R` — Cartografía de severidad (solo fusión)

## Requisitos
* **Google Earth Engine:** Cuenta activa en la plataforma para la ejecución de los *scripts* de la carpeta `00_preprocesamiento_GEE`.
* **Lenguaje R (>= 4.3.0):** Para la ejecución del modelado estadístico y cartográfico.
* **Librerías principales de R:**
  * Tratamiento de datos espaciales y ráster: `terra`, `sf`
  * Modelización de *Machine Learning*: `ranger`, `randomForest`, `Boruta`
  * Manipulación de datos: `dplyr`

## Nota sobre las rutas
Los scripts emplean rutas absolutas locales (por ejemplo, `G:/TFM/...`) 
correspondientes al entorno de trabajo original. Para reproducir los análisis, 
es necesario adaptar estas rutas a la ubicación de los datos en cada equipo.

## Autoría
Trabajo Fin de Máster — Máster Universitario en Tecnologías de la Información 
Geográfica, Universidad de Extremadura.

David San Martín Aguado
