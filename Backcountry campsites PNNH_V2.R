# ANALYSIS MANUSCRIPT BACKCOUNTRY CAMPSITES

#Version 2

# Limpieza del enviroment
ls ()
rm (list = ls())
ls ()


#Librerias
library(patchwork)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(MASS)
library(fitdistrplus)
library(ggeffects)
library(glmmTMB)
library(TMB)
library(car)
library(effects)
library(carData)
library(paletteer)
library(colorBlindness)
library(dplyr)
library(DHARMa)
library(ggthemes)
library(lme4)
library(rstatix)
library(mgcv)
library(ordinal)

#-------------------------------------------------------------------------------------------------------------------------
#-----------cargar datos-----------------

setwd("G:/My Drive/4. 5 lagunas")


datos_distancia <- read.xlsx("datos_campamentos_paisaje.xlsx",   sheet = 1)
datos_usos <- read.xlsx("datos_campamentos_paisaje.xlsx",   sheet = 3)
datos_matriz <- read.xlsx("datos_campamentos_paisaje.xlsx",   sheet = 4)


setwd("G:/My Drive/4. 5 lagunas/Capitulo travesias")
datos_encuesta<- read.xlsx("encuesta_5lag.xlsx",   sheet = 1) 
datos_vegetacion<- read.xlsx("datos_campamentostravesias.xlsx",   sheet = 1)
df_final<- read.xlsx("datos_registro_5lag.xlsx",   sheet = 1)
datos_paisaje<- read.csv("G:/My Drive/4. 5 lagunas/Capitulo travesias/poligonos_metricas.csv", 
                         sep = ",", dec = ".", header = TRUE)



#-------------------------------------------------------------------------------------------------------------------------
#---DISTANCIAS CAMPAMENTOS-------------------


library(sf)
library(purrr)

#shapes
puntos_uso <- st_read("G:/My Drive/4. 5 lagunas/2025_datos sig/curados/usos_camp_subset2.shp")
senda <- st_read("G:/My Drive/4. 5 lagunas/2025_datos sig/curados/sendas_pnnh_travesias_relevadas_3.shp")
hidrografia <- st_read("G:/My Drive/4. 5 lagunas/SIG base/hidrografia_pnnh_lines.shp")

# Estandarizar CRS a metros
senda <- st_transform(senda, st_crs(puntos_uso))
hidrografia <- st_transform(hidrografia, st_crs(puntos_uso))

#NND
calc_micro_dist <- function(p_origen, p_destino, tipo_nombre) {
  if (nrow(p_destino) == 0) return(rep(NA_real_, nrow(p_origen)))
  d_mat <- st_distance(p_origen, p_destino)
  sapply(1:nrow(p_origen), function(i) {
    dists <- as.numeric(d_mat[i, ])
    if (p_origen$name_2[i] == tipo_nombre) dists <- dists[dists > 0.05]
    if (length(dists) == 0) return(NA_real_)
    min(dists, na.rm = TRUE)
  })
}

#-------------------------------------------------------------------------------------------------------------------------
#1A---Calcular distancias a micro escala del campamento


df_puntos_final <- puntos_uso %>%
  group_by(id_camp) %>%
  group_modify(~ {
    p_camp  <- .x %>% filter(name_2 == "camp")
    p_fire  <- .x %>% filter(name_2 == "firepits")
    p_trash <- .x %>% filter(name_2 == "littertrash")
    p_bath  <- .x %>% filter(name_2 == "bath")
    
    .x %>% mutate(
      nnd_camp  = calc_micro_dist(.x, p_camp, "camp"),
      nnd_fire  = calc_micro_dist(.x, p_fire, "firepits"),
      nnd_trash = calc_micro_dist(.x, p_trash, "littertrash"),
      nnd_bath  = calc_micro_dist(.x, p_bath, "bath")
    )
  }) %>%
  ungroup() %>%
  mutate(
    dist_senda = as.numeric(st_distance(geometry, senda[st_nearest_feature(geometry, senda), ], by_element = TRUE)),
    dist_agua  = as.numeric(st_distance(geometry, hidrografia[st_nearest_feature(geometry, hidrografia), ], by_element = TRUE))
  )


#renombrar micro ambiente
df_puntos_final <- df_puntos_final %>%
  rename(ambiente_micro = ambiente)

#amb macro
df_puntos_final <- df_puntos_final %>%
  left_join(
    datos_paisaje %>% select(id, ambiente_macro = ambiente), 
    by = c("id_camp" = "id")
  )

#amb micro
df_puntos_final <- df_puntos_final %>%
  mutate(
    ambiente_micro = case_when(
      ambiente_micro == 1 ~ "Bosque",
      ambiente_micro == 2 ~ "Mallín",
      ambiente_micro == 3 ~ "Roca",
      TRUE ~ as.character(ambiente_micro)
    )
  )



#-------------------------------------------------------------------------------------------------------------------------
#1B---Calcular distancias a macroescala escala del campamento


#Limpiar y asegurar que datos_paisaje sea numérico
datos_paisaje_clean <- datos_paisaje %>%
  mutate(across(where(is.character), ~ gsub(",", ".", .))) %>% 
  mutate(across(c(id, area, LSI,starts_with("elev"), starts_with("slope"), starts_with("orient")), as.numeric))

#ANN
df_ann <- df_puntos_final %>%
  st_drop_geometry() %>%
  group_by(id_camp) %>%
  summarise(
    n_puntos = n(),
    d_obs = mean(pmin(nnd_camp, nnd_fire, nnd_trash, nnd_bath, na.rm = TRUE), na.rm = TRUE)
  ) %>%
  left_join(datos_paisaje_clean %>% select(id, area_campo = area, ambiente, uso), by = c("id_camp" = "id")) %>%
  filter(!is.na(area_campo), area_campo > 0, n_puntos > 2) %>%
  mutate(
    d_exp = 0.5 / sqrt(n_puntos / area_campo),
    ANN = d_obs / d_exp,
    patron = case_when(
      ANN < 1 ~ "Agregado", 
      ANN > 1 ~ "Disperso", 
      TRUE ~ "Aleatorio"
    )
  )

#CÁLCULO DE DISTANCIAS MEDIAS (Internas y Externas)
df_internas_externas <- df_puntos_final %>%
  st_drop_geometry() %>%
  group_by(id_camp) %>%
  summarise(
    # Distancias internas (entre puntos del mismo tipo)
    dist_med_carpas  = mean(nnd_camp[name_2 == "camp"], na.rm = TRUE),
    dist_med_fogones = mean(nnd_fire[name_2 == "firepits"], na.rm = TRUE),
    dist_med_banos   = mean(nnd_bath[name_2 == "bath"], na.rm = TRUE),
    
    # MODIFICACIÓN: Distancias externas calculadas SOLO basándose en las carpas
    dist_med_senda   = mean(dist_senda[name_2 == "camp"], na.rm = TRUE),
    dist_med_agua    = mean(dist_agua[name_2 == "camp"], na.rm = TRUE)
  )

# unir df
df_camp_final <- df_ann %>%
  left_join(df_internas_externas, by = "id_camp") %>%
  left_join(datos_paisaje_clean, by = c("id_camp" = "id")) %>%
  rename(ambiente = ambiente.x) %>%
  select(-matches("\\.y$"))


print(colnames(df_camp_final))






#-------------------------------------------------------------------------------------------------------------------------
#----graficos 


df_grafico_micro <- df_puntos_final %>%
  filter(name_2 %in% c("camp", "firepits", "littertrash", "bath")) %>%
  mutate(
    tipo_impacto = case_when(
      name_2 == "camp" ~ "Carpa",
      name_2 == "firepits" ~ "Fogón",
      name_2 == "littertrash" ~ "Residuos",
      name_2 == "bath" ~ "Baño Informal"
    ),
    ambiente_micro = factor(ambiente_micro, levels = c("Bosque", "Mallín", "Roca"))
  ) %>%
  filter(dist_senda < 100) 

ggplot(df_grafico_micro, aes(x = ambiente_micro, y = dist_senda, fill = tipo_impacto)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              alpha = 0.3, size = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = " ",
    subtitle = " ",
    x = "Ambiente",
    y = "Distancia a la Senda",
    fill = "Usos"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    panel.grid.major.x = element_blank()
  )

df_grafico_completo <- df_puntos_final %>%
  filter(name_2 %in% c("camp", "firepits", "littertrash", "bath")) %>%
  # Pasamos las dos distancias a una sola columna
  pivot_longer(cols = c(dist_senda, dist_agua), 
               names_to = "recurso", 
               values_to = "distancia_metros") %>%
  mutate(
    tipo_impacto = case_when(
      name_2 == "camp" ~ "Carpa",
      name_2 == "firepits" ~ "Fogón",
      name_2 == "littertrash" ~ "Residuos",
      name_2 == "bath" ~ "Baño Informal"
    ),
    recurso = ifelse(recurso == "dist_senda", "Senda Principal", "Fuente de Agua"),
    ambiente_micro = factor(ambiente_micro, levels = c("Bosque", "Mallín", "Roca")),
    ambiente_macro = factor(ambiente_macro, levels = c("bosque", "mallin"), labels = c("Bosque", "Mallín"))
  ) %>%
  filter(distancia_metros <= 150)


ggplot(df_grafico_completo, aes(x = tipo_impacto, y = distancia_metros, fill = recurso)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2, size = 1) +
  scale_fill_manual(values = c("Senda Principal" = "#D95F02", "Fuente de Agua" = "#1B9E77")) +
  labs(x = "Tipo de Impacto", y = "Distancia (metros)", fill = " ") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(df_grafico_completo, aes(x = ambiente_micro, y = distancia_metros, fill = tipo_impacto)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~recurso) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2, size = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = " ", y = "Distancia (metros)", fill = "Uso") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df_grafico_completo, aes(x = tipo_impacto, y = distancia_metros, fill = recurso)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ambiente_macro) + 
  scale_fill_manual(values = c("Senda Principal" = "#D95F02", "Fuente de Agua" = "#1B9E77")) +
  labs(x = "Uso", y = "Distancia (metros)", fill = " ") +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))


#-------------------------------------------------------------------------------------------------------------------------
#PREGUNTA 1
#--------------------------------#
#1A: MACROESCALA


df_macro_modelos <- df_camp_final %>%
  mutate(ambiente = as.factor(ambiente)) %>%
  filter(!is.na(LSI), !is.na(ANN), !is.na(dist_med_carpas))%>%
  filter(dist_med_carpas < 40)


df_macro_final_clean <- df_macro_modelos %>%
  mutate(
    # PASO CLAVE: Recodificar 'muy alto' como 'alto'
    uso_corregido = ifelse(uso.x == "muy_alto", "alto", uso.x),
    
    # Convertir a factor ordenado con los 3 niveles
    intensidad = factor(uso_corregido, levels = c("bajo", "medio", "alto")),
    
    # Asegurar que ambiente sea factor
    ambiente = as.factor(ambiente)
  ) %>%
  # Eliminar filas donde no hay dato de intensidad
  filter(!is.na(intensidad))


# 1. Creamos la variable basada en tu lista
df_macro_final_clean <- df_macro_final_clean %>%
  mutate(
    # Limpiamos el nombre para que no fallen las coincidencias
    nombre_clean = str_to_lower(str_trim(nombre)),
    
    tipo_campamento2 = case_when(
      nombre_clean %in% c("lagunanegra", "jakob", "frey", "ilon1","ilon2", "ilon3") ~ "conservicio",
      nombre_clean %in% c("cab1", "cab2", "mdblanco", "md2", "creton", "creton2", "creton3",
                          "navidadp", "rucaco 1", "rucaco 2", "rucaco 3", "rucaco 4", 
                          "campanile 1", "campanile 2", "llum 1", "llum 2", "navidad2", 
                          "baileywillis", "la chata", "md1", "ricardo", "manolo", "azul", "goye", "jujuy") ~ "sinservicio",
      TRUE ~ "otro" # Por si algún nombre no coincide
    ),
    
    tipo_campamento2 = factor(tipo_campamento2, levels = c("conservicio", "sinservicio")),
    ambiente = as.factor(ambiente),
    intensidad = factor(intensidad, levels = c("bajo", "medio", "alto"))
  )


table(df_macro_final_clean$tipo_campamento2)

table(df_macro_final_clean$tipo_campamento2,
      df_macro_final_clean$ambiente)




#cargamos variables topograficas


library(stringr)


slope_tab <- st_read("G:/My Drive/4. 5 lagunas/2025_datos sig/slope_estadistica.shp") %>% 
  st_drop_geometry() %>%
  dplyr::select(id, slp_mean, slp_stdev)

tri_tab <- st_read("G:/My Drive/4. 5 lagunas/2025_datos sig/tri_estadistica.shp") %>% 
  st_drop_geometry() %>%
  dplyr::select(id, tri_mean, tri_stdev)

# Unir id_camp = id
df_macro_final_clean <- df_macro_final_clean %>%
  left_join(slope_tab, by = c("id_camp" = "id")) %>%
  left_join(tri_tab, by = c("id_camp" = "id"))

# Lista
nombres_lista <- df_macro_final_clean %>% 
  pull(nombre) %>% 
  unique() %>% 
  sort()

print(nombres_lista)


df_macro_final_clean <- df_macro_final_clean %>%
  mutate(zona_campamento = case_when(
    str_detect(str_to_lower(nombre), "creton")    ~ "laguna_creton",
    str_detect(str_to_lower(nombre), "rucaco")    ~ "valle_rucaco",
    str_detect(str_to_lower(nombre), "campanile") ~ "valle_campanile",
    str_detect(str_to_lower(nombre), "llum")      ~ "laguna_llum",
    str_detect(str_to_lower(nombre), "ilon")      ~ "laguna_ilon",
    str_detect(str_to_lower(nombre), "md")        ~ "mallin_matedulce",
    str_detect(str_to_lower(nombre), "cab")       ~ "laguna_cab",
    str_detect(str_to_lower(nombre), "navidad")   ~ "valle_navidad",
    TRUE ~ nombre 
  ))

table(df_macro_final_clean$zona_campamento,df_macro_final_clean$ambiente)

# 1. Definimos la variable respuesta en el df de 30 parches
df_macro_final_clean <- df_macro_final_clean %>%
  group_by(zona_campamento) %>%
  mutate(es_disperso = ifelse(n() > 1, 1, 0)) %>% # 1 si la zona tiene >1 parche
  ungroup()


# Aseguramos que la unidad de agrupamiento sea la zona geográfica
df_macro_final_clean <- df_macro_final_clean %>%
  mutate(
    zona_campamento = as.factor(zona_campamento),
    tipo_campamento2 = factor(tipo_campamento2, levels = c("sinservicio", "conservicio")),
    intensidad = factor(uso_corregido, levels = c("bajo", "medio", "alto"))
  )



#Modelo paisaje 



mod_paisaje <- glmmTMB(area ~ ambiente + dist_med_senda + dist_med_agua, 
                      data = df_macro_final_clean, 
                      family = Gamma(link = "log"),
                      dispformula = ~ ambiente) 

plot(simulateResiduals(mod_paisaje))
summary(mod_paisaje)



mod_paisaje_lsi <- glmmTMB(LSI ~ ambiente + dist_med_senda + dist_med_agua, 
                           data = df_macro_final_clean, 
                           family =  Gamma(link = "log"))

summary(mod_paisaje_lsi)
plot(simulateResiduals(mod_paisaje_lsi))


mod_paisaje_area_2 <- glmmTMB(area ~ ambiente + slp_mean + tri_mean, 
                              data = df_macro_final_clean, 
                              family =  Gamma(link = "log"))

summary(mod_paisaje_area_2)
plot(simulateResiduals(mod_paisaje_area_2))




mod_paisaje_lsi_2 <- glmmTMB(LSI ~ ambiente + slp_mean + tri_mean, 
                             data = df_macro_final_clean, 
                             family =  Gamma(link = "log"))


summary(mod_paisaje_lsi_2)
plot(simulateResiduals(mod_paisaje_lsi_2))



#Modelo de uso 

mod_uso <- glmmTMB(area ~ ambiente + intensidad + tipo_campamento2, 
                   data = df_macro_final_clean, 
                   family =  Gamma(link = "log"))

summary(mod_uso)
plot(simulateResiduals(mod_uso))


mod_uso_lsi <- glmmTMB(LSI ~ ambiente + intensidad + tipo_campamento2, 
                       data = df_macro_final_clean, 
                       family =  Gamma(link = "log"))

summary(mod_uso_lsi)
plot(simulateResiduals(mod_uso_lsi))



#MODELOS DE GESTION

mod_area_tipologia1 <- glmmTMB(area ~ tipo_campamento2 + ambiente, 
                               data = df_macro_final_clean, 
                               family = Gamma(link = "log"),
                               dispformula = ~ tipo_campamento2)


summary(mod_area_tipologia1)
plot(simulateResiduals(mod_area_tipologia1))



# Modelo Gamma con link log (Ideal para áreas)
mod_area_tipologia4 <- glmmTMB(LSI ~ tipo_campamento2 + ambiente, 
                               data = df_macro_final_clean, 
                               family =  Gamma(link = "log"),
                               dispformula = ~ tipo_campamento2)

summary(mod_area_tipologia4)
plot(simulateResiduals(mod_area_tipologia4))



#graficar la salida de los 6 modelos 



library(patchwork)


# PREDICCIONES 

# --- TOPOGRAFÍA (Pendiente) ---
pr_area_topo <- predict_response(mod_paisaje_area_2, terms = c("slp_mean [all]", "ambiente"))
pr_lsi_topo  <- predict_response(mod_paisaje_lsi_2, terms = c("slp_mean [all]", "ambiente"))

# --- INTENSIDAD (Carga Social) ---
pr_area_uso  <- predict_response(mod_uso, terms = c("intensidad", "ambiente"))
pr_lsi_uso   <- predict_response(mod_uso_lsi, terms = c("intensidad", "ambiente"))

# --- TIPOLOGÍA (Gestión) ---
pr_area_typ  <- predict_response(mod_area_tipologia1, terms = c("tipo_campamento2", "ambiente"))
pr_lsi_typ   <- predict_response(mod_area_tipologia4, terms = c("tipo_campamento2", "ambiente"))



preparar_labels <- function(df) {
  df <- as.data.frame(df)
  df %>% mutate(
    group = factor(group, levels = c("bosque", "mallin"), labels = c("Forest", "Meadow")),
    x = case_when(
      x == "bajo" ~ "Low", x == "medio" ~ "Medium", x == "alto" ~ "High",
      x == "conservicio" ~ "With Services", x == "sinservicio" ~ "Without Services",
      TRUE ~ as.character(x)
    )
  )
}

df1 <- preparar_labels(pr_area_topo); df2 <- preparar_labels(pr_area_uso)
df3 <- preparar_labels(pr_area_typ);  df4 <- preparar_labels(pr_lsi_topo)
df5 <- preparar_labels(pr_lsi_uso);   df6 <- preparar_labels(pr_lsi_typ)

# Asegurar orden de factores
df2$x <- factor(df2$x, levels = c("Low", "Medium", "High"))
df5$x <- factor(df5$x, levels = c("Low", "Medium", "High"))
df3$x <- factor(df3$x, levels = c("With Services", "Without Services"))
df6$x <- factor(df6$x, levels = c("With Services", "Without Services"))


paleta <- c("Forest" = "#D95F02", "Meadow" = "#1B9E77")

# --- FILA 1: ÁREA (Magnitud) ---
# Pendiente (Línea)
g1 <- ggplot(df1, aes(x = as.numeric(x), y = predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 1) + scale_color_manual(values = paleta) + scale_fill_manual(values = paleta) +
  coord_cartesian(ylim = c(0, 20000)) + labs(title = " ", y = "Area (m²)", x = " ") + theme_minimal() + theme(legend.position = "none")

# Intensidad (Puntos)
g2 <- ggplot(df2, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.4), size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.4)) +
  scale_color_manual(values = paleta) + coord_cartesian(ylim = c(0, 20000)) + 
  labs(title = " ", y = "", x = " ") + theme_minimal() + theme(legend.position = "none")

# Tipología (Puntos)
g3 <- ggplot(df3, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.4), size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.4)) +
  scale_color_manual(values = paleta) + coord_cartesian(ylim = c(0, 20000)) + 
  labs(title = " ", y = "", x = " ") + theme_minimal() + theme(legend.position = "none")

# --- FILA 2: LSI (Morfología) ---
g4 <- ggplot(df4, aes(x = as.numeric(x), y = predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, color = NA) +
  geom_line(size = 1) + scale_color_manual(values = paleta) + scale_fill_manual(values = paleta) +
  coord_cartesian(ylim = c(1, 2.5)) + labs(y = "LSI Index", x = "Slope (°)") + theme_minimal() + theme(legend.position = "none")

g5 <- ggplot(df5, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.4), size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.4)) +
  scale_color_manual(values = paleta) + coord_cartesian(ylim = c(1, 2.5)) + 
  labs(y = "", x = "Intensity") + theme_minimal() + theme(legend.position = "none")

g6 <- ggplot(df6, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.4), size = 3) + geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.4)) +
  scale_color_manual(values = paleta) + coord_cartesian(ylim = c(1, 2.5)) + 
  labs(y = "", x = "Services") + theme_minimal() + theme(legend.position = "none")


paleta <- c("Forest" = "#D95F02", "Meadow" = "#1B9E77")


# --- FILA 1: ÁREA (Magnitud) ---

# G1: Topografía (Línea) - Ponemos la significancia del Ambiente arriba
g1 <- g1 + annotate("text", x = 12, y = 18000, label = "Env*** / Slope .", size = 3, fontface = "bold") +
  labs(subtitle = NULL) # Quitamos el subtítulo anterior

# G2: Presión Social (Puntos) - Asteriscos sobre el nivel 'High'
g2 <- g2 + annotate("text", x = 2.5, y = 18000, label = "Env*** / High use***", size = 3, fontface = "bold") +
  labs(subtitle = NULL)

# G3: Gestión (Puntos) - Solo significancia de Ambiente
g3 <- g3 + annotate("text", x = 1.6, y = 18000, label = "Env***", size = 3, fontface = "bold") +
  labs(subtitle = NULL)

# --- FILA 2: LSI (Morfología) ---

# G4: Topografía (Línea) - Significancia de Ambiente y Pendiente
g4 <- g4 + annotate("text", x = 12, y = 2.4, label = "Env*** / Slope*", size = 3, fontface = "bold") +
  labs(subtitle = NULL)

# G5: Presión Social (Puntos) - Solo ambiente
g5 <- g5 + annotate("text", x = 2, y = 2.4, label = "Env*", size = 3, fontface = "bold") +
  labs(subtitle = NULL)

# G6: Gestión (Puntos) - Asteriscos sobre la diferencia de Servicios
g6 <- g6 + annotate("text", x = 1, y = 2.3, label = "Serv**", size = 3, fontface = "bold") +
  labs(subtitle = NULL)



mosaico_final <- (g1 | g2 | g3) / (g4 | g5 | g6)

figura_lista <- mosaico_final + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = " ",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
  )

figura_lista











#colinealidad de las variables


tabla_uso_servicios <- table(df_macro_final_clean$intensidad, df_macro_final_clean$tipo_campamento2)
print(tabla_uso_servicios)



library(DescTools)


CramerV(df_macro_final_clean$intensidad, df_macro_final_clean$tipo_campamento2)


# < 0.3: Asociación débil
# 0.3 - 0.5: Asociación moderada 
# > 0.6: Asociación fuerte 


library(car)

mod_vif_test <- lm(as.numeric(area) ~ intensidad + tipo_campamento2 + ambiente, 
                   data = df_macro_final_clean)

vif(mod_vif_test)


# Si GVIF^(1/(2*Df)) es menor a 2.0 (o el GVIF es menor a 5), NO hay colinealidad.


library(corrplot)


df_geografico_corr <- df_macro_final_clean %>%
  st_drop_geometry() %>%
  dplyr::select(slp_mean, tri_mean, dist_med_senda, dist_med_agua) %>%
  filter(complete.cases(.))


matriz_geo <- cor(df_geografico_corr, method = "spearman")


corrplot(matriz_geo, method = "number", type = "upper", 
         tl.col = "black", title = "Correlación Variables Geográficas",
         mar = c(0,0,1,0))

# Si algún valor es > 0.7 o < -0.7, hay una correlación fuerte. 




mod_vif_geo <- lm(as.numeric(area) ~ slp_mean + tri_mean + dist_med_senda + dist_med_agua + ambiente, 
                  data = df_macro_final_clean)


vif_resultados <- vif(mod_vif_geo)
print(vif_resultados)



#-------------------------------------------------------------------------------------------------------------------------

#1B: MICROESCALA

#DISTANCIA AL VECINO MAS CERCANO ENTRE CARPAS

df_hacinamiento_max <- df_puntos_final %>%
  filter(name_2 == "camp") %>%
  filter(nnd_camp >= 0.5 & nnd_camp <= 35) %>% 
  mutate(
    ambiente_micro = factor(ambiente_micro, levels = c("Bosque", "Mallín", "Roca")),
    id_camp = factor(id_camp)
  )

mod_hacinamiento_senda <- gam(nnd_camp ~ ambiente_macro + s(dist_senda, k=3) + s(id_camp, bs="re"), 
                              data = df_hacinamiento_max, 
                              family = scat(link="log"))

summary(mod_hacinamiento_senda)
plot(simulateResiduals(mod_hacinamiento_senda))


plot(predict_response(mod_hacinamiento_senda, terms = c("dist_senda [all]"))) +
  labs(y = "Separación entre carpas (m)", x = "Distancia a la Senda (m)")



pred_nnd_global <- predict_response(mod_hacinamiento_senda, 
                                    terms = "dist_senda [all]")

df_plot_global <- as.data.frame(pred_nnd_global)

ggplot(df_plot_global, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  geom_line(size = 1.2, color = "steelblue") +
  coord_cartesian(ylim = c(0, 30)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  labs(
    title = " ",
    subtitle = " ",
    x = "Distance to trailhead (m)",
    y = "Nearest Neighbor Distance \n between tent sites (m)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )



#DISTANCIA A LA SENDA 


df_puntos_clean <- df_puntos_final %>%
  filter(dist_senda >= 0.05 & dist_senda <= 150) %>%
  filter(dist_agua >= 0.05 & dist_agua <= 150) %>%
  mutate(
    name_2 = factor(name_2),
    ambiente_micro = factor(ambiente_micro),
    ambiente_macro = factor(ambiente_macro),
    id_camp = factor(id_camp)
  )


# 1. Aditivo Macro (Tweedie + Dispersión)
mod_senda_macro_add <- glmmTMB(dist_senda ~ name_2 + ambiente_macro + (1 | id_camp), 
                               data = df_puntos_clean, 
                               family = tweedie(link = "log"),
                               dispformula = ~ name_2)

# 2. Interacción Macro
mod_senda_macro_int <- glmmTMB(dist_senda ~ name_2 * ambiente_macro + (1 | id_camp), 
                               data = df_puntos_clean, 
                               family = tweedie(link = "log"),
                               dispformula = ~ name_2)

# 3. Aditivo Micro
mod_senda_micro_add <- glmmTMB(dist_senda ~ name_2 + ambiente_micro + (1 | id_camp), 
                               data = df_puntos_clean, 
                               family = tweedie(link = "log"),
                               dispformula = ~ name_2)

# 4. Interacción Micro
mod_senda_micro_int <- glmmTMB(dist_senda ~ name_2 * ambiente_micro + (1 | id_camp), 
                               data = df_puntos_clean, 
                               family = tweedie(link = "log"),
                               dispformula = ~ name_2)

modelo_senda<- glmmTMB(dist_senda ~ name_2 + ambiente_macro + (1 | id_camp), 
                       data = df_puntos_clean, 
                       family = Gamma(link = "log"))


plot(simulateResiduals(mod_senda_macro_add))
plot(simulateResiduals(mod_senda_macro_int))
plot(simulateResiduals(mod_senda_micro_add))
plot(simulateResiduals(mod_senda_micro_int))
plot(simulateResiduals(modelo_senda))




summary(mod_senda_macro_add)
summary(mod_senda_macro_int)
summary(mod_senda_micro_add)
summary(mod_senda_micro_int)
summary(modelo_senda)

aic_senda <- AIC(mod_senda_macro_add, mod_senda_macro_int, mod_senda_micro_add, mod_senda_micro_int)
print(aic_senda[order(aic_senda$AIC), ])



#SALIDA DIST SENDA ADD


pred_senda <- predict_response(mod_senda_macro_add, terms = c("name_2", "ambiente_macro"))


df_p1 <- as.data.frame(pred_senda) %>%
  mutate(x = case_when(
    x == "camp" ~ "Campsite",
    x == "firepits" ~ "Firepit",
    x == "bath" ~ "Informal Toilet",
    x == "littertrash" ~ "Trash",
    TRUE ~ as.character(x)
  ),
  group = factor(group, levels = c("bosque", "mallin"), labels = c("Forest", "Meadow")))

ggplot(df_p1, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(0.5), size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.5), size = 1) +
  scale_color_manual(values = c("Forest" = "#D95F02", "Meadow" = "#1B9E77")) +
  coord_cartesian(ylim = c(0, 80)) +
  labs(title = " ", 
       x = "Use type", y = "Distance between use and trailhead (m)", color = " ") +
  theme_minimal() + theme(legend.position = "bottom", text = element_text(size = 11))


#-------------------------------------------------------------------------------------------------------------------------
#PREGUNTA 2: VEGETACION
#-------------------------------------------------------
#unir df

veg_camp_completo<-st_read("G:/My Drive/4. 5 lagunas/2025_datos sig/curados/veg_camp_completo.shp")

library(sf)
library(dplyr)

veg_atrib_preparado <- datos_vegetacion %>%
  mutate(
    across(c(arbusto, subarbusto, herbgram, herb, musgoliquen, sd, sinevidencia), as.numeric),
    cobertura_basal = arbusto + subarbusto + herbgram + herb + musgoliquen,
    # % de Daño
    pct_danio = 100 - sinevidencia
  )

#
veg_spat_completo <- veg_camp_completo %>%
  left_join(veg_atrib_preparado, by = "codigo") %>%
  st_transform(st_crs(senda))



veg_spat_completo <- veg_spat_completo %>%
  mutate(
    dist_senda_real = as.numeric(st_distance(
      geometry, 
      senda[st_nearest_feature(geometry, senda), ], 
      by_element = TRUE
    ))
  )



veg_spat_final <- veg_spat_completo %>%
  group_by(id) %>% # 
  group_modify(~ {
    puntos_este_camp <- puntos_uso %>% filter(id_camp == .y$id)
    
    if (nrow(puntos_este_camp) > 0) {
      idx <- st_nearest_feature(.x, puntos_este_camp)
      .x %>% mutate(
        dist_uso_cercano = as.numeric(st_distance(.x, puntos_este_camp[idx, ], by_element = TRUE)),
        tipo_uso_cercano = puntos_este_camp$name_2[idx]
      )
    } else {
      .x %>% mutate(dist_uso_cercano = NA_real_, tipo_uso_cercano = NA_character_)
    }
  }) %>%
  ungroup()


df_vegetacion_final <- veg_spat_final %>%
  left_join(
    df_camp_final %>% select(id_camp, LSI, ambiente_macro = ambiente),
    by = c("id" = "id_camp")
  ) %>%
  mutate(
    ambiente_macro = as.factor(ambiente_macro),
    pisoteo = as.factor(pisoteo),
    tipo_uso_cercano = as.factor(tipo_uso_cercano)
  )


df_mod_veg <- df_vegetacion_final %>%
  filter(id != 13) %>%
  mutate(
    cobertura_basal = ifelse(cobertura_basal <= 0, 0.1, cobertura_basal),
    pct_danio = ifelse(pct_danio <= 0, 0.1, pct_danio),
    sd = ifelse(sd <= 0, 0.1, sd),
    pisoteo_ord = factor(pisoteo, ordered = TRUE, levels = c("1", "2", "3", "4")),
    id = as.factor(id),
    ambiente_macro = as.factor(ambiente_macro)
  )



n <- nrow(df_vegetacion_final)
df_veg_adj <- df_vegetacion_final %>%
  filter(id != 13) %>%
  mutate(
    # Transformación para que 0 y 1 sean aceptados por la familia Beta
    # y = (y * (n-1) + 0.5) / n
    prop_cobertura = ( (cobertura_basal/100) * (n - 1) + 0.5 ) / n,
    prop_danio = ( (pct_danio/100) * (n - 1) + 0.5 ) / n,
    prop_sd = ( (sd/100) * (n - 1) + 0.5 ) / n, 
    pisoteo_ord = factor(pisoteo, ordered = TRUE, levels = c("1", "2", "3", "4")),
    id = as.factor(id),
    ambiente_macro = as.factor(ambiente_macro)
  )

# correlacion pisoteo cob basal 


cor_test <- cor.test(as.numeric(df_veg_adj$pisoteo), 
                     df_veg_adj$cobertura_basal, 
                     method = "spearman")
print(cor_test)



#Vegetacion modelos 


corr_test <- cor.test(df_veg_adj$distancia, df_veg_adj$dist_uso_cercano, method = "spearman")
print(corr_test)

ggplot(df_veg_adj, aes(x = distancia, y = dist_uso_cercano, color = ambiente_macro)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  labs(title = "Relación entre Gradiente Senda vs. Gradiente Uso",
       x = "Distancia a la Senda (m)", y = "Distancia al Uso Cercano (m)") +
  theme_minimal()


df_veg_adj <- df_veg_adj %>%
  mutate(pisoteo_ord = factor(pisoteo, levels = c("1", "2", "3", "4"), ordered = TRUE),
         id = as.factor(id))



df_veg_adj <- df_veg_adj %>% filter( dist_uso_cercano >= 0)



#write.csv(df_veg_adj, "datos4.csv", row.names = FALSE)
#MODELO PISOTEO 


modelo_pisoteo <- clmm(
  pisoteo_ord ~ distancia * ambiente_macro + (1 | id),
  data = df_veg_adj,
  link = "logit" )


summary(modelo_pisoteo)
sim_res <- simulateResiduals(fittedModel = modelo_pisoteo)
plot(sim_res)

pred_pisoteo <- ggpredict(modelo_pisoteo, terms = c("distancia [all]", "ambiente_macro"))
df_plot <- as.data.frame(pred_pisoteo)



df_plot <- df_plot %>%
  mutate(response.level = case_when(
    response.level == "1" ~ "1. Nulo",
    response.level == "2" ~ "2. Bajo",
    response.level == "3" ~ "3. Medio",
    response.level == "4" ~ "4. Alto"
  ))

ggplot(df_plot, aes(x = x, y = predicted, color = response.level, fill = response.level)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~group, labeller = labeller(group = c("bosque" = "Bosque", 
                                                   "mallin" = "Mallin"))) +
  scale_color_manual(values = c("1. Nulo" = "#2E8B57", "2. Bajo" = "#9ACD32", 
                                "3. Medio" = "#FFA500", "4. Alto" = "#B22222")) +
  scale_fill_manual(values = c("1. Nulo" = "#2E8B57", "2. Bajo" = "#9ACD32", 
                               "3. Medio" = "#FFA500", "4. Alto" = "#B22222")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = " ",
    subtitle = " ",
    x = "Distancia desde la Senda (m)",
    y = "Probabilidad Predicha (%)",
    color = "Nivel de Pisoteo",
    fill = "Nivel de Pisoteo"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )



library(emmeans)


pendientes_pisoteo <- emtrends(modelo_pisoteo, ~ ambiente_macro, var = "distancia")
summary(pendientes_pisoteo)

# Bosque: -0.13 (Pendiente fuerte)
# Mallín: -0.04 (Pendiente suave)

mod_piso_nulo <- clmm(pisoteo_ord ~ 1 + (1 | id), data = df_veg_adj, link = "logit")

# R-cuadrado aproximado (McFadden)
L_null <- -257.72 # logLik de un modelo sin predictores 
L_full <- -219.45 # logLik del summary
R2_McFadden <- 1 - (L_full / L_null)
R2_McFadden


#PRUEBAS CON OTROS MODELOS 

# ¿El pisoteo depende de la proximidad a la carpa/fogon y cambia según el ambiente?
mod_piso_uso <- clmm(pisoteo_ord ~ dist_uso_cercano * ambiente_macro + (1 | id), 
                     data = df_veg_adj, link = "logit")

summary(mod_piso_uso)

# ¿Cuál de los dos gradientes tiene más peso cuando los ponemos juntos?
mod_piso_global <- clmm(pisoteo_ord ~ distancia + dist_uso_cercano + ambiente_macro + (1 | id), 
                        data = df_veg_adj, link = "logit")

summary(mod_piso_global)



# COMPARACION CON EL NULO
# Interpretación: Si Pr(>Chisq) es < 0.05, tus predictores son significativos.

print(anova(mod_piso_nulo, modelo_pisoteo))

print(anova(mod_piso_nulo, mod_piso_uso))

print(anova(mod_piso_nulo, mod_piso_global))







#-------------------------------------------------------------------------------------------------------------------------
#PREGUNTA 3: PERCEPCIONES
#-------------------------------------------------------


datos_encuesta <- datos_encuesta %>%
  mutate(
    ciudad_limpia = str_trim(tolower(procedencia_ciudad)),
    procedencia = case_when(
      ciudad_limpia %in% c("bariloche") ~ "local",
      TRUE ~ "no-local"
    ),
    procedencia = as.factor(procedencia)
  )


ggplot(datos_encuesta, aes(x = fct_infreq(procedencia_provincia), fill = procedencia)) +
  geom_bar() +
  coord_flip() + 
  theme_minimal() +
  labs(
    title = "Procedencia de los visitantes por Provincia",
    x = "Provincia",
    y = "Cantidad de Encuestados",
    fill = "Categoría"
  )


datos_encuesta$desdecuando_encuestado <- factor(
  datos_encuesta$desdecuando_encuestado, 
  levels = c(
    "Recientemente", 
    "En los últimos 2 años", 
    "En los últimos 5 años", 
    "Desde siempre"
  ),
  ordered = TRUE
)


datos_encuesta$frecuencia_visita_encuestado <- factor(
  datos_encuesta$frecuencia_visita_encuestado,
  levels = c(
    "Es la primera vez", 
    "Cada tanto", 
    "Una vez por año", 
    "Más de una vez por año"
  ),
  ordered = TRUE
)

levels(datos_encuesta$desdecuando_encuestado)
levels(datos_encuesta$frecuencia_visita_encuestado)
datos_encuesta$af_firepits <- as.factor(datos_encuesta$af_firepits)


#3.0 PATRONES DE USO 
#------------------------------------------------------------------------------------------------------------------------


datos_acampe_largo <- datos_encuesta %>%
  pivot_longer(
    cols = starts_with("acampe_dia"), 
    names_to = "dia_columna",         
    values_to = "lugar",              
    values_drop_na = TRUE             
  ) %>%
  mutate(
    lugar = str_trim(lugar),
    lugar_limpio = case_when(
      lugar %in% c("No acampamos", "no acampamos", "Ninguno", "Hotel", "Vivac en el filo del Cerro Capitan") ~ "Otro",
      TRUE ~ lugar
    ),
  
    dia_label = str_replace(dia_columna, "acampe_dia", "Día ")
  ) %>%
  filter(!is.na(lugar_limpio), lugar_limpio != "")


unique(datos_acampe_largo$lugar_limpio)


ggplot(datos_acampe_largo, aes(x = dia_label, fill = lugar_limpio)) +
  geom_bar(position = "stack", color = "white", size = 0.1) +
  labs(title = "Campamentos elegidos según el día de travesía",
       subtitle = "Categorías 'Hotel', 'Vivac' y 'No acampamos' agrupadas en 'Otro'",
       x = "Día",
       y = "Cantidad de personas",
       fill = "Sitio de acampe") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 0)
  )


frecuencia_sitios <- datos_acampe_largo %>%
  count(lugar_limpio) %>%
  arrange(desc(n))

ggplot(frecuencia_sitios, aes(x = reorder(lugar_limpio, n), y = n, fill = n)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  labs(title = "Presión de uso total por sitio",
       x = "Sitio de acampe",
       y = "Total de pernoctes registrados") +
  theme_minimal() +
  theme(legend.position = "none")




sitios_principales <- frecuencia_sitios %>%
  filter(lugar_limpio != "Otro") %>%
  filter(n >= n[lugar_limpio == "Arroyo La Chata"])


nombres_principales <- sitios_principales$lugar_limpio


datos_temporal <- datos_acampe_largo %>%
  filter(lugar_limpio %in% nombres_principales) %>%
  count(dia_label, lugar_limpio) %>%
  mutate(dia_label = factor(dia_label, levels = c("Día 1", "Día 2", "Día 3", "Día 4", "Día 5", "Día 6")))



ggplot(datos_temporal, aes(x = dia_label, y = n, color = lugar_limpio, group = lugar_limpio)) +
  geom_line(size = 1.2, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Uso de sitios principales a lo largo de la travesía",
       subtitle = "Cantidad de pernoctes por día en los campamentos más frecuentados",
       x = "Día de la travesía",
       y = "Cantidad de visitantes",
       color = "Sitio de acampe") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  ) +
  facet_wrap(~lugar_limpio, scales = "free_y")



datos_agrupados_temporal <- datos_acampe_largo %>%
  mutate(lugar_zona = case_when(
    str_detect(lugar_limpio, "Laguna CAB") ~ "Laguna CAB",
    str_detect(lugar_limpio, "Cretón") ~ "Cretón",
    TRUE ~ lugar_limpio
  )) %>%
  filter(lugar_zona %in% c("Laguna CAB", "Cretón", 
                           "Laguna Ilón", "Laguna Negra", "Arroyo La Chata", 
                           "Mallín del Mate Dulce/ Mallín de Las Vueltas")) %>%
  count(dia_label, lugar_zona) %>%
  mutate(dia_label = factor(dia_label, levels = c("Día 1", "Día 2", "Día 3", "Día 4", "Día 5", "Día 6")))


head(datos_agrupados_temporal)

ggplot(datos_agrupados_temporal, aes(x = dia_label, y = n, color = lugar_zona, group = lugar_zona)) +
  geom_line(size = 1.3, alpha = 0.8) +
  geom_point(size = 3.5) +
  scale_color_brewer(palette = "Set1") + 
  labs(title = " ",
       subtitle = " ",
       x = " ",
       y = "Cantidad de pernoctes",
       color = "Zona de acampe") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10)
  )



datos_procedencia <- datos_acampe_largo %>%
  filter(lugar_limpio != "Otro") %>%
  count(lugar_limpio, procedencia) %>%
  filter(!is.na(procedencia))


ggplot(datos_procedencia, aes(x = reorder(lugar_limpio, n), y = n, fill = procedencia)) +
  geom_bar(stat = "identity", position = "fill") + # "fill" hace que todas las barras midan 1 (100%)
  coord_flip() +
  scale_y_continuous(labels = scales::percent) + # Eje en porcentaje
  scale_fill_manual(values = c("local" = "#2E8B57", "no-local" = "#DAA520")) +
  labs(title = " ",
       subtitle = " ",
       x = "Sitio de acampe",
       y = "Porcentaje de visitantes",
       fill = "Procedencia") +
  theme_minimal()




#--------------------------------------

#armar grafico de costo distancia-altimetria 


library(tidyverse)
library(openxlsx)

# 1. Leer el archivo de costo
setwd("G:/My Drive/4. 5 lagunas")
datos_costo <- read.xlsx("datos_campamentos_paisaje.xlsx", sheet = 5)

# --- CORRECCIÓN DEL ERROR ---
# Desduplicamos los nombres de las columnas del Excel (el segundo 'nombre' pasa a ser 'nombre.1')
names(datos_costo) <- make.unique(names(datos_costo))

# 2. Unir con los datos de acampe
datos_acampe_unido <- datos_acampe_largo %>%
  mutate(nombre_costo = case_when(
    lugar_limpio == "Rancho Manolo"                                ~ "manolo",
    lugar_limpio == "Laguna Negra"                                 ~ "lagunanegra",
    lugar_limpio == "Arroyo La Chata"                              ~ "la chata",
    lugar_limpio == "Laguna CAB (llegada)"                         ~ "cab1",
    lugar_limpio == "Laguna CAB (después de rodearla)"             ~ "cab2",
    lugar_limpio == "Mallín del Mate Dulce/ Mallín de Las Vueltas" ~ "md1",
    lugar_limpio == "Laguna Cretón"                                ~ "creton",
    lugar_limpio == "Pozones Cretón"                               ~ "creton2",
    lugar_limpio == "Laguna Jujuy"                                 ~ "jujuy",
    lugar_limpio == "Mallin de Ricardo"                            ~ "ricardo",
    lugar_limpio == "Laguna Ilón"                                  ~ "ilon2",
    TRUE ~ NA_character_
  )) %>%
  left_join(
    datos_costo %>% select(nombre, any_of(c("cost", "distancia_km", "elev_gain", "elevacion_ganada_m"))),
    by = c("nombre_costo" = "nombre")
  )

# 3. Resumen por campamento
camps_resumen <- datos_acampe_unido %>%
  filter(!is.na(cost)) %>%
  group_by(lugar_limpio, cost) %>%
  summarise(total_pernoctes = n(), .groups = "drop") %>%
  filter(total_pernoctes > 0) %>%
  arrange(cost) %>%
  mutate(
    distancia_km = cost / 1000,
    id_orden     = row_number(),
    referencia   = paste0(id_orden, "=", lugar_limpio)
  )

# 4. Curva suavizada (spline)
curva_suave <- as.data.frame(
  spline(camps_resumen$distancia_km, camps_resumen$total_pernoctes, n = 300)
) %>%
  mutate(y = pmax(0, y)) # Evita valores negativos

# 5. Definición de zonas por día
zonas_dias <- tibble(
  etiqueta_dia = c("Day 1", "Day 2", "Day 3", "Day 4 & 5"),
  xmin         = c(5.0,  13.0, 22.0, 29.0),
  xmax         = c(13.0, 22.0, 29.0, 35.5),
  x_centro     = (xmin + xmax) / 2,
  fill_color   = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
)

leyenda_camps <- paste(camps_resumen$referencia, collapse = "  |  ")

# 6. Gráfico Final
ggplot() +
  # Franjas de los días
  geom_rect(data = zonas_dias, 
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = etiqueta_dia),
            alpha = 0.10, show.legend = FALSE) +
  
  # Textos "Day 1", "Day 2", etc.
  geom_text(data = zonas_dias, 
            aes(x = x_centro, y = 100, label = etiqueta_dia, color = etiqueta_dia),
            fontface = "bold", size = 6, lineheight = 0.9, show.legend = FALSE) +
  
  # Curva spline
  geom_line(data = curva_suave, aes(x = x, y = y), 
            color = "grey40", linewidth = 1, linetype = "dashed") +
  
  # Puntos negros numerados
  geom_point(data = camps_resumen, aes(x = distancia_km, y = total_pernoctes),
             color = "#1E2A38", size = 10) +
  
  geom_text(data = camps_resumen, aes(x = distancia_km, y = total_pernoctes, label = id_orden),
            color = "white", fontface = "bold", size = 5) +
  
  # Escalas y estética
  scale_fill_manual(values = zonas_dias$fill_color) +
  scale_color_manual(values = zonas_dias$fill_color) +
  scale_x_continuous(
    breaks = seq(5, 35, by = 5),
    labels = function(x) paste0(x, " km"),
    limits = c(5, 35.5)
  ) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    x = "Trail distance (km)",
    y = "Total number of recorded overnight stays",
    caption = paste0("Campsites: ", leyenda_camps)
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14),
    plot.caption = element_text(size = 9, color = "grey35", hjust = 0, margin = margin(t = 13))
  )







#Sin suavizar 




# Gráfico con línea directa entre los datos crudos (sin suavizado)
ggplot() +
  # 1. Franjas de fondo de los días
  geom_rect(data = zonas_dias, 
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = etiqueta_dia),
            alpha = 0.10, show.legend = FALSE) +
  
  # 2. Textos "Day 1", "Day 2", etc.
  geom_text(data = zonas_dias, 
            aes(x = x_centro, y = 100, label = etiqueta_dia, color = etiqueta_dia),
            fontface = "bold", size = 6, lineheight = 0.9, show.legend = FALSE) +
  
  # 3. Línea DIRECTA conectando los datos crudos reales (sin suavizado)
  geom_line(data = camps_resumen, 
            aes(x = distancia_km, y = total_pernoctes), 
            color = "grey40", linewidth = 0.9, linetype = "dashed") +
  
  # 4. Puntos oscuros con el número de campamento dentro
  geom_point(data = camps_resumen, 
             aes(x = distancia_km, y = total_pernoctes),
             color = "#1E2A38", size = 10) +
  
  geom_text(data = camps_resumen, 
            aes(x = distancia_km, y = total_pernoctes, label = id_orden),
            color = "white", fontface = "bold", size = 5) +
  
  # 5. Escalas y leyendas
  scale_fill_manual(values = zonas_dias$fill_color) +
  scale_color_manual(values = zonas_dias$fill_color) +
  scale_x_continuous(
    breaks = seq(5, 35, by = 5),
    labels = function(x) paste0(x, " km"),
    limits = c(5, 35.5)
  ) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(
    x = "Trail distance (km)",
    y = "Total number of recorded overnight stays",
    caption = paste0("Campsites: ", leyenda_camps)
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14),
    plot.caption = element_text(size = 9, color = "grey35", hjust = 0, margin = margin(t = 13))
  )




#--------------------------------------------------------------
#3.1 ARMADO DE RUTAS POSIBLES TRAVESIA 5 LAGUNAS 


#tabla de visitantes por campamento


# visitantes por dia por campamento
tabla_visitantes <- datos_acampe_largo %>%
  mutate(
    tamañogrupo = as.numeric(tamañogrupo),
    tamañogrupo = replace_na(tamañogrupo, 1)
  ) %>%
  group_by(campamento = lugar_limpio, dia = dia_label) %>%
  summarise(
    total_visitantes = sum(tamañogrupo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dia, desc(total_visitantes))


print(tabla_visitantes, n = Inf)


# write.csv(tabla_visitantes, "ruta_preferida.csv", row.names = FALSE)

# visitantes totales por campamento 
tabla_total_campamento <- datos_acampe_largo %>%
  mutate(
    tamañogrupo = as.numeric(tamañogrupo),
    tamañogrupo = replace_na(tamañogrupo, 1)
  ) %>%
  group_by(campamento = lugar_limpio) %>%
  summarise(
    total_visitantes = sum(tamañogrupo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_visitantes))

# Ver la tabla completa
print(tabla_total_campamento, n = Inf)


#write.csv(tabla_total_campamento, "ruta_preferida2.csv", row.names = FALSE)





#rutas realizadas por los visitantes para estimar la ruta preferida 

rutas_por_visitante <- datos_encuesta %>%
  mutate(
    id  = row_number(),
    tamañogrupo = as.numeric(tamañogrupo),
    tamañogrupo = replace_na(tamañogrupo, 1)
  ) %>%
  pivot_longer(
    cols = starts_with("acampe_dia"),
    names_to = "dia_col",
    values_to = "lugar",
    values_drop_na = TRUE
  ) %>%
  mutate(
    lugar = str_trim(lugar),
    lugar_limpio = case_when(
      lugar %in% c("No acampamos", "no acampamos", "Ninguno", "Hotel", "Vivac en el filo del Cerro Capitan") ~ "Otro",
      TRUE ~ lugar
    )
  ) %>%
  filter(!is.na(lugar_limpio), lugar_limpio != "", lugar_limpio != "Otro") %>%
  # Unir los campamentos de cada persona con una flecha " -> "
  group_by(id) %>%
  summarise(
    ruta_completa  = paste(lugar_limpio, collapse = " -> "),
    noches_totales = n(),
    permanenciatravesia      = first(permanenciatravesia),
    tamanogrupo    = first(tamañogrupo),
    .groups = "drop"
  )


#write.csv(rutas_por_visitante, "ruta_preferida_visitantes2.csv", row.names = FALSE)


#-----------------------------------------------------------------------------------

############################################
#Tabla resumen topografia rutas



library(tidyverse)
library(openxlsx)

setwd("G:/My Drive/4. 5 lagunas")
archivo_excel <- "datos_campamentos_paisaje.xlsx"

# ==============================================================================
# 1. CARGA DE DATOS FÍSICOS (Hoja 6: Distancias, desniveles y tiempos por nodo)
# ==============================================================================
df_fijos <- read.xlsx(archivo_excel, sheet = 6) %>% rename_with(make.unique)

df_tramos_base <- df_fijos %>%
  mutate(
    distancia_km_earth      = as.numeric(distancia_km_earth),
    elev_ganada             = as.numeric(elev_ganada),
    elev_perdida            = as.numeric(elev_perdida),
    total_tiempo_tramo_min  = as.numeric(total_tiempo_tramo_min),
    id_18                   = row_number()
  ) %>%
  select(id_18, nombre_costo, distancia_km_earth, elev_ganada, elev_perdida, total_tiempo_tramo_min)

# ==============================================================================
# 2. CARGA DE LA TABLA DE RUTAS (Hoja 7)
# ==============================================================================
df_rutas_7 <- read.xlsx(archivo_excel, sheet = 7) %>% rename_with(make.unique)
names(df_rutas_7)[1:5] <- c("ruta", "dia_pernocte", "nombre_costo", "costo_acumulado", "visitantes_acumulados")

# Filtramos las rutas que necesitas
df_rutas_limpias <- df_rutas_7 %>%
  mutate(ruta = str_trim(tolower(ruta))) %>%
  filter(!is.na(ruta), ruta %in% c("optima", "conservadora", "rapida", "ruta1", "ruta2", "ruta3", "parque")) %>%
  mutate(
    dia_pernocte = as.numeric(dia_pernocte),
    ruta_nombre = case_when(
      ruta == "optima"       ~ "1. Ruta Óptima (6-7h)",
      ruta == "conservadora" ~ "2. Ruta Conservadora (4-5h)",
      ruta == "rapida"       ~ "3. Ruta Rápida (7-9h)",
      ruta == "ruta1"        ~ "4. Ruta 5 Días (Clásica)",
      ruta == "ruta2"        ~ "5. Ruta 4 Días (Directa CAB)",
      ruta == "ruta3"        ~ "6. Ruta 3 Días (Rápida Exprés)",
      ruta == "parque"       ~ "7. Ruta Parque Nacional",
      TRUE ~ ruta
    )
  )

# ==============================================================================
# 3. UNIÓN Y CÁLCULO DE MÉTRICAS (TIEMPO SOLO EN MINUTOS)
# ==============================================================================
# Calculamos los acumulados por campamento en la base física
df_tramos_acum <- df_tramos_base %>%
  mutate(
    dist_acum_nodo     = cumsum(distancia_km_earth),
    elev_pos_acum_nodo = cumsum(elev_ganada),
    tiempo_acum_nodo   = cumsum(total_tiempo_tramo_min)
  )

# Cruzamos cada pernocte de cada ruta con sus acumulados
tabla_rutas_minutos <- df_rutas_limpias %>%
  left_join(df_tramos_acum, by = "nombre_costo") %>%
  group_by(ruta_nombre) %>%
  arrange(dia_pernocte) %>%
  mutate(
    # Métricas acumuladas en ese punto de la ruta
    Distancia_Acum_km    = round(dist_acum_nodo, 2),
    Desnivel_Pos_Acum_m  = round(elev_pos_acum_nodo, 1),
    Tiempo_Acum_min      = round(tiempo_acum_nodo, 1),
    
    # Métricas de ESE tramo/jornada específica (restando el día anterior)
    Distancia_Tramo_km   = round(Distancia_Acum_km - lag(Distancia_Acum_km, default = 0), 2),
    Desnivel_Pos_Tramo_m = round(Desnivel_Pos_Acum_m - lag(Desnivel_Pos_Acum_m, default = 0), 1),
    Tiempo_Tramo_min     = round(Tiempo_Acum_min - lag(Tiempo_Acum_min, default = 0), 1)
  ) %>%
  ungroup() %>%
  # Seleccionamos las columnas conservando solo minutos para el tiempo
  select(
    Ruta                 = ruta_nombre,
    Dia                  = dia_pernocte,
    Campamento_Destino   = nombre_costo,
    Distancia_Tramo_km,
    Distancia_Acum_km,
    Desnivel_Pos_Tramo_m,
    Desnivel_Pos_Acum_m,
    Tiempo_Tramo_min,
    Tiempo_Acum_min
  )

# Ver la tabla completa en consola
print(tabla_rutas_minutos, n = Inf)



#write.csv(tabla_rutas_minutos, "tabla_rutas_completa2.csv", row.names = FALSE)


##########################################################


## 4 ESCENARIOS + PERFIL TOPOGRAFICO 


setwd("G:/My Drive/4. 5 lagunas")
archivo_excel <- "datos_campamentos_paisaje.xlsx"

#Perfil Topográfico
df_perfil_relieve <- read.csv("perfil_senda.csv") %>%
  rename(distancia_m = distance, elevacion = elevation) %>%
  mutate(
    distancia_m  = as.numeric(distancia_m),
    elevacion    = as.numeric(elevacion),
    distancia_km = distancia_m / 1000 
  ) %>%
  filter(!is.na(distancia_km), !is.na(elevacion))

max_km_real <- max(df_perfil_relieve$distancia_km)


zonas_dias <- tibble(
  etiqueta_dia = c("Day 1", "Day 2", "Day 3", "Day 4 & 5"),
  xmin         = c(0.0,  13.0, 22.0, 29.0),
  xmax         = c(13.0, 22.0, 29.0, max_km_real),
  x_centro     = (xmin + xmax) / 2,
  fill_color   = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
)


df_fijos <- read.xlsx(archivo_excel, sheet = 6) %>% rename_with(make.unique)

df_18_nodos_v2 <- df_fijos %>%
  mutate(
    distancia_km_earth = as.numeric(distancia_km_earth),
    km_acumulado       = cumsum(distancia_km_earth),
    id_18              = row_number(),
    total_visitantes   = replace_na(as.numeric(total_visitantes), 0),
    referencia         = paste0(id_18, "=", nombre_costo)
  ) %>%
  select(id_18, nombre_costo, km_acumulado, total_visitantes, referencia)

#Pampalinda no es campamento es punto de llegada
df_18_sin_pampa_h6 <- df_18_nodos_v2 %>% filter(id_18 != 18) 
total_global_h6    <- sum(df_18_sin_pampa_h6$total_visitantes, na.rm = TRUE)

df_real_h6_pct <- df_18_sin_pampa_h6 %>%
  mutate(
    porcentaje_visitantes = round((total_visitantes / total_global_h6) * 100, 1)
  )

# Interpolacion cota-campamento
df_real_h6_pct$elev_perfil <- approx(
  x = df_perfil_relieve$distancia_km, 
  y = df_perfil_relieve$elevacion, 
  xout = df_real_h6_pct$km_acumulado
)$y

# Leyenda
fila1_sp <- paste(df_real_h6_pct$referencia[1:6], collapse = "   |   ")
fila2_sp <- paste(df_real_h6_pct$referencia[7:12], collapse = "   |   ")
fila3_sp <- paste(df_real_h6_pct$referencia[13:17], collapse = "   |   ")
leyenda_h6 <- paste("Campamentos:", fila1_sp, fila2_sp, fila3_sp, sep = "\n")


#4 escenarios 

df_rutas_raw <- read.xlsx(archivo_excel, sheet = 7) %>% rename_with(make.unique)
names(df_rutas_raw)[1:5] <- c("ruta", "dia_pernocte", "nombre_costo", "costo_acumulado", "visitantes_acumulados")


orden_4_rutas <- c(
  "3 Days",
  "4 Days (Directa CAB)",
  "4 Days (Salida Cretón)",
  "5 Days"
)

df_nuevas_rutas <- df_rutas_raw %>%
  mutate(ruta_limpia = str_replace_all(str_trim(tolower(ruta)), " ", "")) %>%
  filter(ruta_limpia %in% c("ruta1", "ruta2", "ruta3", "ruta4")) %>%
  mutate(
    visitantes_acumulados = as.numeric(visitantes_acumulados),
    dia_pernocte          = as.numeric(dia_pernocte)
  ) %>%
  left_join(df_18_nodos_v2 %>% select(nombre_costo, km_acumulado), by = "nombre_costo") %>%
  mutate(
    ruta_label = case_when(
      ruta_limpia == "ruta3" ~ "3 Days",
      ruta_limpia == "ruta2" ~ "4 Days (Directa CAB)",
      ruta_limpia == "ruta4" ~ "4 Days (Salida Cretón)",
      ruta_limpia == "ruta1" ~ "5 Days"
    ),
    ruta_label = factor(ruta_label, levels = orden_4_rutas)
  ) %>%
  arrange(ruta_label, dia_pernocte)

# Etapas continuas (Ribbon)
df_tramos_etapas <- df_nuevas_rutas %>%
  group_by(ruta_label) %>%
  mutate(
    km_inicio = lag(km_acumulado, default = 0),
    km_fin    = km_acumulado
  ) %>%
  ungroup()

df_ribbon_etapas <- df_tramos_etapas %>%
  rowwise() %>%
  do({
    r <- .
    df_perfil_relieve %>%
      filter(distancia_km >= r$km_inicio & distancia_km <= r$km_fin) %>%
      mutate(
        ruta_label   = factor(r$ruta_label, levels = orden_4_rutas),
        dia_pernocte = as.factor(r$dia_pernocte)
      )
  }) %>%
  ungroup()


#perfil topografico grafico 

fig_perfil_h6_publicacion <- ggplot() +
  geom_rect(data = zonas_dias, aes(xmin = xmin, xmax = xmax, ymin = 700, ymax = Inf, fill = etiqueta_dia), alpha = 0.08, show.legend = FALSE) +
  geom_ribbon(data = df_perfil_relieve, aes(x = distancia_km, ymin = 750, ymax = elevacion), fill = "#900C3F", alpha = 0.15) +
  geom_line(data = df_perfil_relieve, aes(x = distancia_km, y = elevacion), color = "#900C3F", linewidth = 1.0) +
  geom_line(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil), color = "grey35", linewidth = 0.9, linetype = "dotted") +
  geom_point(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil, size = porcentaje_visitantes), color = "#1E2A38", alpha = 0.9) +
  geom_text(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil, label = id_18), color = "white", fontface = "bold", size = 3.0, show.legend = FALSE) +
  
  scale_fill_manual(values = zonas_dias$fill_color) +
  scale_size_continuous(range = c(2.5, 10), breaks = c(0, 5, 10, 15, 20), name = "Visitors use (%)") +
  scale_x_continuous(breaks = seq(5, 35, by = 5), labels = function(x) paste0(x, " km"), limits = c(0, max_km_real), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(800, 2000, by = 300), limits = c(750, 2250), labels = function(y) paste0(y, " m")) +
  
  labs(
    title = " ", subtitle = " ",
    x = "Trail Network", y = "Elevation (m)",
    caption = leyenda_h6
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    axis.title         = element_text(face = "bold", size = 11),
    axis.text          = element_text(size = 14),
    legend.position    = "right",
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 14),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption       = element_text(size = 10.0, color = "grey35", hjust = 0, lineheight = 1.2, margin = margin(t = 12))
  )

print(fig_perfil_h6_publicacion)


#grafico combinado perfil + escenerarios


p1_preferida <- ggplot() +
  geom_rect(data = zonas_dias, aes(xmin = xmin, xmax = xmax, ymin = 700, ymax = Inf, fill = etiqueta_dia), alpha = 0.08, show.legend = FALSE) +
  geom_ribbon(data = df_perfil_relieve, aes(x = distancia_km, ymin = 750, ymax = elevacion), fill = "#900C3F", alpha = 0.15) +
  geom_line(data = df_perfil_relieve, aes(x = distancia_km, y = elevacion), color = "#900C3F", linewidth = 0.9) +
  geom_line(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil), color = "grey35", linewidth = 0.8, linetype = "dotted") +
  geom_point(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil, size = porcentaje_visitantes), color = "#1E2A38", alpha = 0.9) +
  geom_text(data = df_real_h6_pct, aes(x = km_acumulado, y = elev_perfil, label = id_18), color = "white", fontface = "bold", size = 2.8, show.legend = FALSE) +
  
  scale_fill_manual(values = zonas_dias$fill_color) +
  scale_size_continuous(range = c(2.5, 9), name = "Visitors use (%)") +
  scale_x_continuous(breaks = seq(5, 35, by = 5), limits = c(0, max_km_real), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(800, 2000, by = 400), limits = c(750, 2250), labels = function(y) paste0(y, " m")) +
  
  labs(title = " ", x = NULL, y = " ") +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey92"),
    axis.title         = element_text(face = "bold", size = 11),
    axis.text          = element_text(size = 10),
    axis.text.x        = element_blank(), 
    legend.position    = "right",
    plot.title         = element_text(face = "bold", size = 12)
  )


p_escenarios <- ggplot() +
  geom_rect(data = zonas_dias, aes(xmin = xmin, xmax = xmax, ymin = 700, ymax = Inf, fill = etiqueta_dia), alpha = 0.08, show.legend = FALSE) +
  geom_ribbon(data = df_perfil_relieve, aes(x = distancia_km, ymin = 750, ymax = elevacion), fill = "grey85", alpha = 0.4) +
  geom_ribbon(data = df_ribbon_etapas, aes(x = distancia_km, ymin = 750, ymax = elevacion, fill = dia_pernocte), alpha = 0.55) +
  geom_line(data = df_perfil_relieve, aes(x = distancia_km, y = elevacion), color = "#900C3F", linewidth = 0.8) +
  
  
  geom_point(data = df_nuevas_rutas %>% filter(nombre_costo != "pampalinda"), 
             aes(x = km_acumulado, y = approx(df_perfil_relieve$distancia_km, df_perfil_relieve$elevacion, xout = (df_nuevas_rutas %>% filter(nombre_costo != "pampalinda"))$km_acumulado)$y), 
             color = "black", size = 3) +
  
  
  facet_wrap(~ruta_label, ncol = 1) +
  
  scale_fill_brewer(palette = "Pastel1", name = "Day") +
  scale_x_continuous(breaks = seq(5, 35, by = 5), labels = function(x) paste0(x), limits = c(0, max_km_real), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(800, 2000, by = 400), limits = c(750, 2250), labels = function(y) paste0(y, " m")) +
  
  labs(x = "Trail Network (km)", y = "Elevation (m)") +
  theme_minimal() +
  theme(
    strip.text         = element_text(face = "bold", size = 11),
    strip.background   = element_rect(fill = "grey95", color = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "grey92"),
    axis.title         = element_text(face = "bold", size = 11),
    axis.text          = element_text(size = 14),
    legend.position    = "right"
  )

#unir
figura_integrada_5_paneles <- p1_preferida / p_escenarios + 
  plot_layout(heights = c(1.0, 4.0)) + 
  plot_annotation(
    title = " ",
    caption = leyenda_h6,
    theme = theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0),
      plot.caption = element_text(size = 12.0, color = "grey35", hjust = 0, lineheight = 1.2, margin = margin(t = 10))
    )
  )

print(figura_integrada_5_paneles)







#---------------------------------

#datos registro de trekking PNNH - 5 LAGUNAS-


setwd("G:/My Drive/4. 5 lagunas/Capitulo travesias")
df_final<- read.xlsx("datos_registro_5lag.xlsx",   sheet = 1)



df_limpio <- df_final %>%
  mutate(
    f_inicio = make_date(year = as.numeric(year_inicio), 
                         month = as.numeric(mes_inicio), 
                         day = as.numeric(dia_inicio)),
    f_final  = make_date(year = as.numeric(year_final), 
                         month = as.numeric(mes_final), 
                         day = as.numeric(dia_final))
  ) %>%
  filter(!is.na(f_inicio), !is.na(f_final)) %>%
  mutate(
    temp_inicio = pmin(f_inicio, f_final),
    temp_final  = pmax(f_inicio, f_final)
  ) %>%
  select(-f_inicio, -f_final) %>%
  rename(f_inicio = temp_inicio, f_final = temp_final) %>%
  mutate(duracion_dias = as.numeric(f_final - f_inicio) + 1) %>%
  filter(duracion_dias >= 1 & duracion_dias <= 12)


if (exists("datos_encuesta")) {
  pernoctes_para_unir <- df_limpio %>%
    filter(!is.na(mail)) %>%
    select(mail, f_inicio, f_final, duracion_dias) %>%
    distinct(mail, .keep_all = TRUE)
  
  datos_encuesta <- datos_encuesta %>%
    mutate(mail = str_to_lower(str_trim(mail))) %>%
    left_join(pernoctes_para_unir, by = "mail")
}


df_ocupacion_diaria <- df_limpio %>%
  rowwise() %>%
  mutate(dia = list(seq(f_inicio, f_final, by = "day"))) %>%
  unnest(dia) %>%
  ungroup() %>%
  group_by(dia) %>%
  summarise(
    personas_simultaneas = sum(visitantes, na.rm = TRUE), 
    n_grupos = n(),
    .groups = "drop"
  ) %>%
  complete(dia = seq(min(dia), max(dia), by = "day"), 
           fill = list(personas_simultaneas = 0, n_grupos = 0))

df_temporadas <- df_ocupacion_diaria %>%
  filter(month(dia) %in% c(11, 12, 1, 2, 3, 4)) %>%
  mutate(
    temporada = if_else(
      month(dia) >= 11,
      paste0(year(dia), "-", year(dia) + 1),
      paste0(year(dia) - 1, "-", year(dia))
    ),
    dia_mes = day(dia),
    mes_nombre = factor(
      month(dia), 
      levels = c(4, 3, 2, 1, 12, 11), 
      labels = c("Abr", "Mar", "Feb", "Ene", "Dic", "Nov")
    ),
    dia_mes_ficticio = as.Date(paste0(
      if_else(month(dia) >= 11, "2000-", "2001-"),
      sprintf("%02d-%02d", month(dia), day(dia))
    ))
  ) %>%
  filter(temporada %in% c("2020-2021", "2021-2022", "2022-2023", "2023-2024", "2024-2025"))


g1_heatmap <- ggplot(df_temporadas, aes(x = dia_mes, y = mes_nombre, fill = personas_simultaneas)) +
  geom_tile(color = "white", linewidth = 0.15) + 
  facet_wrap(~temporada, ncol = 2) +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#fee08b", "#f46d43", "#d73027", "#67001f"),
    name = "Visitors"
  ) +
  scale_x_continuous(breaks = seq(5, 30, by = 5)) +
  labs(
    title = " ",
    subtitle = " ",
    x = "Day",
    y = "Month"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13)
  )

print(g1_heatmap)




g2_curva <- ggplot(df_temporadas, aes(x = dia_mes_ficticio, y = personas_simultaneas)) +
  geom_area(fill = "#d95f02", alpha = 0.45) +
  geom_line(color = "#d95f02", linewidth = 0.8) +
  facet_wrap(~temporada, ncol = 1, scales = "free_y") +
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b",
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = " ",
    subtitle = " ",
    x = "Season",
    y = "Visitors"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 11, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 13)
  )

print(g2_curva)




# dias con mas de 75 visitantes 



dias_mas_75 <- df_temporadas %>%
  filter(personas_simultaneas >= 75 & personas_simultaneas <= 150)


# ¿Cuántos días por temporada se superaron los 100 visitantes?
resumen_75 <- dias_mas_75 %>%
  group_by(temporada) %>%
  summarise(
    dias_con_mas_75 = n(),
    .groups = "drop"
  )

# ¿En qué meses ocurrió?
meses_75 <- dias_mas_75 %>%
  group_by(temporada, mes_nombre) %>%
  summarise(
    cantidad_dias = n(),
    .groups = "drop"
  )

print(resumen_75)

print(meses_75)






#-------------------------------------------------------------------------------------------------------------------------
#3A correlacion 



#SUBSET DE LA ENCUESTA


usos_conteo <- datos_usos %>%
  separate_rows(name_2, sep = "/") %>% 
  group_by(id_camp, name_2) %>%
  summarise(cantidad = n(), .groups = "drop") %>%
  pivot_wider(names_from = name_2, values_from = cantidad, values_fill = 0) %>%
  rename_with(~paste0("n_", .), -id_camp)


referencia_sig_completa <- df_camp_final %>%
  left_join(usos_conteo, by = "id_camp")

cat("Referencia SIG completa creada con métricas de paisaje y conteos de uso.\n")



traductor_nombres <- data.frame(
  nombre_encuesta = c("Laguna Negra", "Laguna CAB (llegada)", "Laguna CAB (llegada) ", 
                      "Laguna CAB (después de rodearla)", "Arroyo La Chata",
                      "Mallín del Mate Dulce/ Mallín de Las Vueltas", "Pozones Cretón",
                      "Laguna Cretón", "Laguna Ilón", "Mallin Goye", "Laguna Jujuy", 
                      "Laguna Azul", "Rancho Manolo", "Mallin de Ricardo"),
  id_sig = c(19, 4, 4, 5, 3, 8, 9, 10, 16, 35, 12, 20, 18, 14) 
)


df_encuesta_espacial <- datos_encuesta %>%
  pivot_longer(
    cols = starts_with("acampe_dia"), 
    names_to = "nro_noche", 
    values_to = "nombre_encuesta"
  ) %>%
  filter(!is.na(nombre_encuesta), nombre_encuesta != "") %>%
  left_join(traductor_nombres, by = "nombre_encuesta") %>%
  left_join(referencia_sig_completa, by = c("id_sig" = "id_camp"))


df_tesis_final <- df_encuesta_espacial %>%
  mutate(across(starts_with("sig_cant_"), ~replace_na(., 0)))

colnames(df_tesis_final)


df_tesis_final<- df_tesis_final %>%
  mutate(tipo_visitante = case_when(
    frecuencia_visita_encuestado == "Es la primera vez" ~ "nuevo",
    TRUE ~ "repite"
  )) %>%
  mutate(tipo_visitante = as.factor(tipo_visitante))


#write.csv(df_tesis_final, "datos7.csv", row.names = FALSE)

#CORRELACION
 

pisoteo_por_sitio <- veg_spat_final %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(
    sig_pisoteo_medio_sitio = mean(as.numeric(as.character(pisoteo)), na.rm = TRUE),
    .groups = "drop"
  )

df_exposicion_media <- df_tesis_final %>%
  st_drop_geometry() %>%
  left_join(pisoteo_por_sitio, by = c("id_sig" = "id")) %>%
  group_by(id) %>% 
  summarise(
    p_fire  = first(as.numeric(af_firepits)),
    p_trash = first(as.numeric(af_trash)),
    p_bath  = first(as.numeric(af_informalbaths)),
    p_crowd = first(as.numeric(af_crowding)),
    p_veg   = first(as.numeric(af_trampling)), 
    sig_n_fire  = mean(n_firepits, na.rm = TRUE),
    sig_d_fire  = mean(dist_med_fogones, na.rm = TRUE),
    sig_n_trash = mean(n_littertrash, na.rm = TRUE),
    sig_n_bath  = mean(n_bath, na.rm = TRUE),
    sig_d_bath  = mean(dist_med_banos, na.rm = TRUE),
    sig_n_camp  = mean(n_camp, na.rm = TRUE),
    sig_d_camp  = mean(dist_med_carpas, na.rm = TRUE),
    sig_pisoteo = mean(sig_pisoteo_medio_sitio, na.rm = TRUE),
    
    .groups = "drop"
  ) 


run_spearman <- function(var_p, var_sig, etiqueta) {
  datos_par <- df_exposicion_media %>% 
    filter(!is.na(!!sym(var_p)), !is.na(!!sym(var_sig)))
  
  if(nrow(datos_par) > 5) {
    test <- cor.test(datos_par[[var_p]], datos_par[[var_sig]], method = "spearman")
    return(data.frame(Impacto = etiqueta, 
                      Rho = round(test$estimate, 3), 
                      P_Value = round(test$p.value, 3),
                      N_Efectivo = nrow(datos_par)))
  } else {
    return(data.frame(Impacto = etiqueta, Rho = NA, P_Value = NA, N_Efectivo = nrow(datos_par)))
  }
}


tabla_final_158 <- bind_rows(
  run_spearman("p_fire",  "sig_n_fire",  "Fogones (Cantidad)"),
  run_spearman("p_fire",  "sig_d_fire",  "Fogones (Distancia)"),
  run_spearman("p_trash", "sig_n_trash", "Basura (Cantidad)"),
  run_spearman("p_bath",  "sig_n_bath",  "Baños (Cantidad)"),
  run_spearman("p_bath",  "sig_d_bath",  "Baños (Distancia)"),
  run_spearman("p_crowd", "sig_n_camp",  "Hacinamiento (Cantidad)"),
  run_spearman("p_crowd", "sig_d_camp",  "Hacinamiento (Distancia)"),
  run_spearman("p_veg",   "sig_pisoteo", "Vegetación (Nivel Pisoteo)")
)


print(tabla_final_158)




# en vez del promedio usamos MIN Y MAX


df_exposicion_maxima <- df_tesis_final %>%
  st_drop_geometry() %>%
  left_join(pisoteo_por_sitio, by = c("id_sig" = "id")) %>%
  group_by(id) %>% 
  summarise(
    p_fire = first(as.numeric(af_firepits)),
    p_trash = first(as.numeric(af_trash)),
    p_bath = first(as.numeric(af_informalbaths)),
    p_crowd = first(as.numeric(af_crowding)),
    p_veg = first(as.numeric(af_trampling)),
    
    #MAX
    sig_n_fire = max(n_firepits, na.rm = TRUE),
    sig_n_trash = max(n_littertrash, na.rm = TRUE),
    sig_n_bath = max(n_bath, na.rm = TRUE),
    sig_n_camp = max(n_camp, na.rm = TRUE),
    sig_pisoteo = max(sig_pisoteo_medio_sitio, na.rm = TRUE),
    
    #MIN
    sig_d_fire = min(dist_med_fogones, na.rm = TRUE),
    sig_d_bath = min(dist_med_banos, na.rm = TRUE),
    sig_d_camp = min(dist_med_carpas, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("sig_"), ~ifelse(is.infinite(.), NA, .)))


run_spearman_max <- function(var_p, var_sig, etiqueta) {
  datos_par <- df_exposicion_maxima %>% 
    filter(!is.na(!!sym(var_p)), !is.na(!!sym(var_sig)))
  
  if(nrow(datos_par) > 5) {
    test <- cor.test(datos_par[[var_p]], datos_par[[var_sig]], method = "spearman")
    return(data.frame(Impacto = etiqueta, 
                      Rho_Max = round(test$estimate, 3), 
                      P_Value_Max = round(test$p.value, 3),
                      N_Efectivo = nrow(datos_par)))
  } else {
    return(data.frame(Impacto = etiqueta, Rho_Max = NA, P_Value_Max = NA, N_Efectivo = nrow(datos_par)))
  }
}


tabla_final_maxima <- bind_rows(
  run_spearman_max("p_fire",  "sig_n_fire",  "Fogones (Cantidad)"),
  run_spearman_max("p_fire",  "sig_d_fire",  "Fogones (Distancia)"),
  run_spearman_max("p_trash", "sig_n_trash", "Basura (Cantidad)"),
  run_spearman_max("p_bath",  "sig_n_bath",  "Baños (Cantidad)"),
  run_spearman_max("p_bath",  "sig_d_bath",  "Baños (Distancia)"),
  run_spearman_max("p_crowd", "sig_n_camp",  "Hacinamiento (Cantidad)"),
  run_spearman_max("p_crowd", "sig_d_camp",  "Hacinamiento (Distancia)"),
  run_spearman_max("p_veg",   "sig_pisoteo", "Vegetación (Nivel Pisoteo)")
)

print(tabla_final_maxima)





#-----correcccion min y max, sin el promedio x campamento SE PIERDE EL EFECTO DE LA SIGNIFICANCIA. MEJOR USAR 1 ANTERIOR




df_internas_externas <- df_puntos_final %>%
  st_drop_geometry() %>%
  group_by(id_camp) %>%
  summarise(
    # Promedios (para el análisis de medias)
    dist_med_carpas  = mean(nnd_camp[name_2 == "camp"], na.rm = TRUE),
    dist_med_fogones = mean(nnd_fire[name_2 == "firepits"], na.rm = TRUE),
    dist_med_banos   = mean(nnd_bath[name_2 == "bath"], na.rm = TRUE),
    dist_med_senda   = mean(dist_senda[name_2 == "camp"], na.rm = TRUE),
    dist_med_agua    = mean(dist_agua[name_2 == "camp"], na.rm = TRUE),
    
    # VERDADEROS MÍNIMOS REALES (El punto más cercano dentro de cada campamento)
    min_d_carpas     = suppressWarnings(min(nnd_camp[name_2 == "camp"], na.rm = TRUE)),
    min_d_fogones    = suppressWarnings(min(nnd_fire[name_2 == "firepits"], na.rm = TRUE)),
    min_d_banos      = suppressWarnings(min(nnd_bath[name_2 == "bath"], na.rm = TRUE)),
    
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("min_d_"), ~ ifelse(is.infinite(.), NA, .)))


df_camp_final <- df_ann %>%
  left_join(df_internas_externas, by = "id_camp") %>%
  left_join(datos_paisaje_clean, by = c("id_camp" = "id")) %>%
  rename(ambiente = ambiente.x) %>%
  select(-matches("\\.y$"))



usos_conteo <- datos_usos %>%
  separate_rows(name_2, sep = "/") %>% 
  group_by(id_camp, name_2) %>%
  summarise(cantidad = n(), .groups = "drop") %>%
  pivot_wider(names_from = name_2, values_from = cantidad, values_fill = 0) %>%
  rename_with(~paste0("n_", .), -id_camp)


referencia_sig_completa <- df_camp_final %>%
  left_join(usos_conteo, by = "id_camp")


traductor_nombres <- data.frame(
  nombre_encuesta = c("Laguna Negra", "Laguna CAB (llegada)", "Laguna CAB (llegada) ", 
                      "Laguna CAB (después de rodearla)", "Arroyo La Chata",
                      "Mallín del Mate Dulce/ Mallín de Las Vueltas", "Pozones Cretón",
                      "Laguna Cretón", "Laguna Ilón", "Mallin Goye", "Laguna Jujuy", 
                      "Laguna Azul", "Rancho Manolo", "Mallin de Ricardo"),
  id_sig = c(19, 4, 4, 5, 3, 8, 9, 10, 16, 35, 12, 20, 18, 14) 
)


df_encuesta_espacial <- datos_encuesta %>%
  pivot_longer(
    cols = starts_with("acampe_dia"), 
    names_to = "nro_noche", 
    values_to = "nombre_encuesta"
  ) %>%
  filter(!is.na(nombre_encuesta), nombre_encuesta != "") %>%
  left_join(traductor_nombres, by = "nombre_encuesta") %>%
  left_join(referencia_sig_completa, by = c("id_sig" = "id_camp"))

df_tesis_final <- df_encuesta_espacial %>%
  mutate(across(starts_with("sig_cant_"), ~replace_na(., 0))) %>%
  mutate(tipo_visitante = as.factor(case_when(
    frecuencia_visita_encuestado == "Es la primera vez" ~ "nuevo",
    TRUE ~ "repite"
  )))



# Pisoteo
pisoteo_por_sitio <- veg_spat_final %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(
    sig_pisoteo_medio_sitio = mean(as.numeric(as.character(pisoteo)), na.rm = TRUE),
    .groups = "drop"
  )

# EXPOSICIÓN MÁXIMA / MÍNIMA REAL
df_exposicion_maxima <- suppressWarnings(
  df_tesis_final %>%
    st_drop_geometry() %>%
    left_join(pisoteo_por_sitio, by = c("id_sig" = "id")) %>%
    group_by(id) %>% 
    summarise(
      # Percepción
      p_fire  = first(as.numeric(af_firepits)),
      p_trash = first(as.numeric(af_trash)),
      p_bath  = first(as.numeric(af_informalbaths)),
      p_crowd = first(as.numeric(af_crowding)),
      p_veg   = first(as.numeric(af_trampling)),
      
      # VERDADEROS MÁXIMOS DE CANTIDAD DEL VIAJE
      sig_n_fire  = max(n_firepits, na.rm = TRUE),
      sig_n_trash = max(n_littertrash, na.rm = TRUE),
      sig_n_bath  = max(n_bath, na.rm = TRUE),
      sig_n_camp  = max(n_camp, na.rm = TRUE),
      sig_pisoteo = max(sig_pisoteo_medio_sitio, na.rm = TRUE),
      
      # VERDADEROS MÍNIMOS DE DISTANCIA DEL VIAJE (MÁXIMA CERCANÍA)
      sig_d_fire  = min(min_d_fogones, na.rm = TRUE),
      sig_d_bath  = min(min_d_banos, na.rm = TRUE),
      sig_d_camp  = min(min_d_carpas, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(across(starts_with("sig_"), ~ ifelse(is.infinite(.), NA, .)))
)

#SPEARMAN
run_spearman_max <- function(var_p, var_sig, etiqueta) {
  datos_par <- df_exposicion_maxima %>% 
    filter(!is.na(!!sym(var_p)), !is.na(!!sym(var_sig)))
  
  if(nrow(datos_par) > 5) {
    test <- cor.test(datos_par[[var_p]], datos_par[[var_sig]], method = "spearman", exact = FALSE)
    return(data.frame(
      Impacto     = etiqueta, 
      Rho_Max     = round(test$estimate, 3), 
      P_Value_Max = round(test$p.value, 3),
      N_Efectivo  = nrow(datos_par)
    ))
  } else {
    return(data.frame(Impacto = etiqueta, Rho_Max = NA, P_Value_Max = NA, N_Efectivo = nrow(datos_par)))
  }
}

# TABLA FINAL
tabla_final_maxima <- bind_rows(
  run_spearman_max("p_fire",  "sig_n_fire",  "Fogones (Cantidad Máxima)"),
  run_spearman_max("p_fire",  "sig_d_fire",  "Fogones (Distancia Mínima)"),
  run_spearman_max("p_trash", "sig_n_trash", "Basura (Cantidad Máxima)"),
  run_spearman_max("p_bath",  "sig_n_bath",  "Baños (Cantidad Máxima)"),
  run_spearman_max("p_bath",  "sig_d_bath",  "Baños (Distancia Mínima)"),
  run_spearman_max("p_crowd", "sig_n_camp",  "Hacinamiento (Cantidad Máxima Carpas)"),
  run_spearman_max("p_crowd", "sig_d_camp",  "Hacinamiento (Distancia Mínima Carpas)"),
  run_spearman_max("p_veg",   "sig_pisoteo", "Vegetación (Máximo Nivel Pisoteo)")
)

print(tabla_final_maxima)




#PRUEBA presencia / ausencia basura en campamentos 




df_trash_unico <- df_tesis_final %>%
  group_by(id, procedencia, tipo_visitante) %>%
  summarise(
    af_trash_m = mean(as.numeric(as.character(af_trash)), na.rm = TRUE),
    
    # ANY es más robusto que max para presencia/ausencia
    # Si al menos un registro tiene basura > 0, devuelve TRUE (1), sino FALSE (0)
    presencia_basura = as.numeric(any(n_littertrash > 0, na.rm = TRUE)),
    .groups = "drop"
  )


cor.test(df_trash_unico$af_trash_m, df_trash_unico$presencia_basura, method = "spearman")






#-------------------------------------------------------------------------------------------------------------------------

#3B FACTORES SOCIALES 

datos_encuesta <- datos_encuesta %>%
  mutate(tipo_visitante = case_when(
    frecuencia_visita_encuestado == "Es la primera vez" ~ "nuevo",
    TRUE ~ "repite"
  ))%>%
  mutate(as.factor(edad_encuestado),
         as.factor(educacion_encuestado),
         as.factor(tipo_visitante),
         as.factor(genero_encuestado), 
         as.factor(desdecuando_encuestado))

library(ordinal)
library(brant)
library(car)
library(DescTools)


#write.csv(datos_encuesta, "datos8.csv", row.names = FALSE)


#MODELOS
m_fire <- polr(as.factor(af_firepits) ~  procedencia + edad_encuestado + educacion_encuestado + tipo_visitante + genero_encuestado, data = datos_encuesta, Hess = TRUE)
m_bath <- polr(as.factor(af_informalbaths) ~ procedencia + edad_encuestado + educacion_encuestado + tipo_visitante + genero_encuestado, data = datos_encuesta, Hess = TRUE)
m_crowd <- polr(as.factor(af_crowding) ~ procedencia + edad_encuestado + educacion_encuestado + tipo_visitante + genero_encuestado, data = datos_encuesta, Hess = TRUE)
m_veg <- polr(as.factor(af_trampling) ~ procedencia + edad_encuestado + educacion_encuestado + tipo_visitante + genero_encuestado, data = datos_encuesta, Hess = TRUE)


#LISTA
modelos_finales <- list(Fogones = m_fire, Baños = m_bath, Hacinamiento = m_crowd, vegetacion =m_veg)

#VALIDACION
for(nombre in names(modelos_finales)) {
  cat("\n\n--- VALIDACIÓN MODELO:", nombre, "---\n")
  mod <- modelos_finales[[nombre]]
  
  # Brant Test
  try(print(brant(mod)), silent = TRUE)
  
  # VIF (proxy con modelo lineal)
  # Usamos la variable respuesta numéricamente solo para el VIF
  y_num <- as.numeric(as.factor(mod$model[,1]))
  print(vif(lm(y_num ~ ., data = mod$model[,-1])))
  
  # McFadden R2
  cat("R2 McFadden:", PseudoR2(mod, which = "McFadden"), "\n")
}


#SUMMARYS
for (nombre in names(modelos_finales)) {
  cat("\n\n============================================\n")
  cat("--- RESUMEN:", nombre, "---\n")
  print(summary(modelos_finales[[nombre]]))
}

# 3. Comparación de ajuste (AIC)
cat("\n\n============================================\n")
cat("--- COMPARACIÓN AIC: ¿N o D? ---\n")
resultados_aic <- data.frame(
  Uso = c("Fogones", "Baños", "Hacinamiento", "vegetacion"),
  AIC = c(AIC(m_fire), AIC(m_bath), AIC(m_crowd), AIC (m_veg))
)


print(resultados_aic)


#Modelos mas simples 


# Ajustar versiones simples
m_fire_4 <- polr(as.factor(af_firepits) ~ procedencia + tipo_visitante + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_fire_1 <- polr(as.factor(af_firepits) ~ procedencia + tipo_visitante, data = datos_encuesta, Hess = TRUE)
m_fire_2 <- polr(as.factor(af_firepits) ~ procedencia + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_fire_3 <- polr(as.factor(af_firepits) ~ procedencia, data = datos_encuesta, Hess = TRUE)


m_bath_4 <- polr(as.factor(af_informalbaths) ~ procedencia + tipo_visitante + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_bath_1 <- polr(as.factor(af_informalbaths) ~ procedencia + tipo_visitante, data = datos_encuesta, Hess = TRUE)
m_bath_2 <- polr(as.factor(af_informalbaths) ~ procedencia + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_bath_3 <- polr(as.factor(af_informalbaths) ~ procedencia, data = datos_encuesta, Hess = TRUE)

m_crow_4 <- polr(as.factor(af_crowding) ~ procedencia + tipo_visitante + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_crow_1 <- polr(as.factor(af_crowding) ~ procedencia + tipo_visitante, data = datos_encuesta, Hess = TRUE)
m_crow_2 <- polr(as.factor(af_crowding) ~ procedencia + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_crow_3 <- polr(as.factor(af_crowding) ~ procedencia, data = datos_encuesta, Hess = TRUE)

m_veg_4 <- polr(as.factor(af_trampling) ~ procedencia + tipo_visitante + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_veg_1 <- polr(as.factor(af_trampling) ~ procedencia + tipo_visitante, data = datos_encuesta, Hess = TRUE)
m_veg_2 <- polr(as.factor(af_trampling) ~ procedencia + edad_encuestado, data = datos_encuesta, Hess = TRUE)
m_veg_3 <- polr(as.factor(af_trampling) ~ procedencia, data = datos_encuesta, Hess = TRUE)



modelos_finales <- list(
  m_fire_1=m_fire_1, m_fire_2=m_fire_2, m_fire_3=m_fire_3, m_fire_4=m_fire_4,
  m_bath_1=m_bath_1, m_bath_2=m_bath_2, m_bath_3=m_bath_3, m_bath_4=m_bath_4,
  m_crow_1=m_crow_1, m_crow_2=m_crow_2, m_crow_3=m_crow_3, m_crow_4=m_crow_4,
  m_veg_1=m_veg_1, m_veg_2=m_veg_2, m_veg_3=m_veg_3, m_veg_4=m_veg_4
)

# summarys
for (nombre in names(modelos_finales)) {
  cat("\n\n======================================================\n")
  cat("--- RESUMEN DEL MODELO:", nombre, "---\n")
  cat("======================================================\n")
  mod <- modelos_finales[[nombre]]
  ctable <- coef(summary(mod))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  tabla_p <- cbind(ctable, "p value" = round(p, 5))
  
  print(tabla_p)
  cat("\nR2 McFadden:", round(PseudoR2(mod, which = "McFadden"), 4), "\n")
}



conteo_edad <- table(datos_encuesta$edad_encuestado)
print("Distribución de la muestra por EDAD:")
print(conteo_edad)


print(table(datos_encuesta$procedencia))



library(MASS)
library(DescTools)
library(dplyr)


percepciones_lista <- list(
  "Fogones" = "af_firepits", 
  "Baños" = "af_informalbaths", 
  "Hacinamiento" = "af_crowding", 
  "Vegetacion" = "af_trampling"
)

resumen_comparativo <- data.frame()

for(nom in names(percepciones_lista)) {
  var_y <- percepciones_lista[[nom]]
  datos_limpios <- library(MASS)
  library(DescTools)
  library(dplyr)
  
  
  percepciones_lista <- list(
    "Fogones" = "af_firepits", 
    "Baños" = "af_informalbaths", 
    "Hacinamiento" = "af_crowding", 
    "Vegetacion" = "af_trampling"
  )
  
  resumen_comparativo <- data.frame()
  
  for(nom in names(percepciones_lista)) {
    var_y <- percepciones_lista[[nom]]
    

    datos_limpios <- datos_encuesta%>%
      filter(!is.na(!!sym(var_y)), !is.na(procedencia), !is.na(tipo_visitante), !is.na(edad_encuestado)) %>%
      mutate(y = as.factor(!!sym(var_y)))
    

    m1 <- polr(y ~ procedencia + tipo_visitante, data = datos_limpios, Hess = TRUE)
    m2 <- polr(y ~ procedencia + edad_encuestado, data = datos_limpios, Hess = TRUE)
    m3 <- polr(y ~ procedencia, data = datos_limpios, Hess = TRUE)
    m4 <- polr(y ~ procedencia + tipo_visitante + edad_encuestado, data = datos_limpios, Hess = TRUE)
    

    resumen_comparativo <- rbind(resumen_comparativo, data.frame(
      Percepcion = nom,
      AIC_M1_Proc_Tipo = AIC(m1),
      AIC_M2_Proc_Edad = AIC(m2),
      AIC_M3_Solo_Proc = AIC(m3),
      AIC_M4_Completo  = AIC(m4)
    ))
  }
  
  print("--- TABLA AIC ---")
  print(resumen_comparativo) %>%
    filter(!is.na(!!sym(var_y)), !is.na(procedencia), !is.na(tipo_visitante), !is.na(edad_encuestado)) %>%
    mutate(y = as.factor(!!sym(var_y)))
  
  # Ajustamos los 4 modelos
  m1 <- polr(y ~ procedencia + tipo_visitante, data = datos_limpios, Hess = TRUE)
  m2 <- polr(y ~ procedencia + edad_encuestado, data = datos_limpios, Hess = TRUE)
  m3 <- polr(y ~ procedencia, data = datos_limpios, Hess = TRUE)
  m4 <- polr(y ~ procedencia + tipo_visitante + edad_encuestado, data = datos_limpios, Hess = TRUE)
  

  resumen_comparativo <- rbind(resumen_comparativo, data.frame(
    Percepcion = nom,
    AIC_M1_Proc_Tipo = AIC(m1),
    AIC_M2_Proc_Edad = AIC(m2),
    AIC_M3_Solo_Proc = AIC(m3),
    AIC_M4_Completo  = AIC(m4)
  ))
}

print("--- TABLA AIC  ---")
print(resumen_comparativo)


library(MASS)
library(dplyr)
library(DescTools)
library(brant)
library(car)

# MODELOS SELECCIONADOS
# Fogones: M4 (Completo) | Baños: M3 (Solo Proc) | Hacinamiento: M3 (Solo Proc) | Vegetación: M2 (Proc+Edad)

modelos_ganadores <- list(
  Fogones = m_fire_4,
  Banos = m_bath_3,
  Hacinamiento = m_crow_3,
  Vegetacion = m_veg_2
)

# 2. VALIDACIÓN FINAL Y LRT
resumen_tesis_percepcion <- data.frame()

for(nombre in names(modelos_ganadores)) {
  cat("\n======================================================")
  cat("\nANALIZANDO GANADOR PARA:", toupper(nombre))
  cat("\n======================================================\n")
  
  mod <- modelos_ganadores[[nombre]]
  
  # A. P-VALORES Y COEFICIENTES
  ctable <- coef(summary(mod))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  cat("\n[COEFICIENTES Y P-VALORES]\n")
  print(round(cbind(ctable, "p value" = p), 5))
  
  # B. TEST DE RAZÓN DE VEROSIMILITUD (LRT vs Nulo)
  # ¿Es el modelo significativamente mejor que un modelo sin variables?

  formula_string <- as.character(formula(mod))
  mod_nulo <- polr(as.formula(paste(formula_string[2], "~ 1")), data = mod$model, Hess = TRUE)
  lrt <- anova(mod_nulo, mod)
  p_lrt <- lrt$`Pr(Chi)`[2]
  

  r2_mcf <- PseudoR2(mod, which = "McFadden")
  

  p_brant <- NA
  if(ncol(mod$model) > 1) {
    b_res <- try(brant(mod), silent = TRUE)
    if(!inherit(b_res, "try-error")) p_brant <- b_res[1,3]
  }
  

  resumen_tesis_percepcion <- rbind(resumen_tesis_percepcion, data.frame(
    Variable = nombre,
    AIC = AIC(mod),
    P_LRT_vs_Nulo = p_lrt,
    R2_McFadden = r2_mcf,
    P_Brant_Omnibus = p_brant
  ))
}

print(resumen_tesis_percepcion)


library(MASS)
library(dplyr)
library(DescTools)
library(brant)


percepciones <- c("af_firepits", "af_informalbaths", "af_crowding", "af_trampling")
nombres_perc <- c("FOGONES", "BAÑOS", "HACINAMIENTO", "VEGETACION")


ranking_aic_final <- data.frame()


for(i in 1:length(percepciones)) {
  v_y <- percepciones[i]
  n_y <- nombres_perc[i]
  
  cat("\n\n################################################################")
  cat("\n  ANALIZANDO PERCEPCIÓN:", n_y)
  cat("\n################################################################\n")
  
  df_clean <- datos_encuesta %>%
    filter(!is.na(!!sym(v_y)), !is.na(procedencia), !is.na(tipo_visitante), !is.na(edad_encuestado)) %>%
    mutate(y = as.factor(!!sym(v_y)))
  
  forms <- list(
    M1 = "y ~ procedencia + tipo_visitante",
    M2 = "y ~ procedencia + edad_encuestado",
    M3 = "y ~ procedencia",
    M4 = "y ~ procedencia + tipo_visitante + edad_encuestado"
  )
  
  # Sub-bucle para correr y validar las 4 versiones
  for(m_name in names(forms)) {
    cat("\n--- MODELO", m_name, "---")
    
    try({
      # A. AJUSTE
      mod <- polr(as.formula(forms[[m_name]]), data = df_clean, Hess = TRUE)
      
      # B. SUMMARY CON P-VALUES
      ctable <- coef(summary(mod))
      p_vals <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
      cat("\n[Coeficientes]\n")
      print(round(cbind(ctable, "p value" = p_vals), 5))
      
      # C. VALIDACIÓN (LRT vs Nulo)
      mod_nulo <- update(mod, . ~ 1)
      lrt <- anova(mod_nulo, mod)
      p_lrt <- lrt$`Pr(Chi)`[2]
      
      # D. R2 y AIC
      r2_mcf <- PseudoR2(mod, which = "McFadden")
      aic_val <- AIC(mod)
      
      cat("\nAIC:", round(aic_val, 2), "| LRT p-val:", round(p_lrt, 5), "| R2 McFadden:", round(r2_mcf, 4), "\n")
      
      # E. BRANT TEST
      cat("[Brant Test Omnibus p-val]:")
      p_brant <- tryCatch({ brant(mod)[1,3] }, error = function(e) NA)
      cat(round(p_brant, 4), "\n")
      
      # Guardar en ranking
      ranking_aic_final <- rbind(ranking_aic_final, data.frame(
        Percepcion = n_y,
        Modelo = m_name,
        AIC = aic_val,
        P_LRT = p_lrt,
        R2 = r2_mcf,
        Brant = p_brant
      ))
      
    }, silent = FALSE)
  }
}


print(ranking_aic_final %>% arrange(Percepcion, AIC))



summary(m_fire_4)

summary(m_bath_1)





#------------------------------------
#PREGUNTA 4: ACEPTABILIDAD
#------------------------------------


library(MASS)
library(lmtest)
library(sandwich)

# tipo_visitante 
df_encuesta_158 <- datos_encuesta %>% 
  mutate(tipo_visitante = case_when(
    frecuencia_visita_encuestado == "Es la primera vez" ~ "nuevo",
    TRUE ~ "repite"
  ),
  tipo_visitante = as.factor(tipo_visitante),
  procedencia = as.factor(procedencia),
  educacion = as.factor(educacion_encuestado),
  edad_encuestado= as.factor(edad_encuestado))


analizar_aceptabilidad_158 <- function(df, var_accion) {
  

  datos <- df %>% 
    filter(!is.na(!!sym(var_accion))) %>%
    mutate(y = as.factor(!!sym(var_accion))) %>%
    na.omit(select(., y, procedencia, tipo_visitante, educacion, edad_encuestado))
  
 
  m <- polr(y ~ procedencia + tipo_visitante + educacion + edad_encuestado, 
            data = datos, Hess = TRUE)
  
  print(summary(m)) # Mostramos summary porque no necesitamos robustez por cluster (ya es 1 fila = 1 persona)
  return(m)
}

modelos_accion_158 <- list(
  cupos = analizar_aceptabilidad_158(df_encuesta_158, "action_cuposvisit"),
  fogones = analizar_aceptabilidad_158(df_encuesta_158, "action_nofire"),
  infra = analizar_aceptabilidad_158(df_encuesta_158, "action_campsite")
)





library(MASS)


analizar_aceptabilidad_158 <- function(df, var_accion) {
  datos <- df %>% 
    filter(!is.na(!!sym(var_accion))) %>%
    mutate(y = as.factor(!!sym(var_accion))) %>%
    select(y, procedencia, tipo_visitante, educacion, edad_encuestado) %>%
    na.omit()
  
  vars_validas <- names(datos)[sapply(datos, function(x) length(unique(x)) > 1)]
  vars_in_model <- setdiff(vars_validas, "y")
  
  formula_mod <- as.formula(paste("y ~", paste(vars_in_model, collapse = "+")))
  
  m <- polr(formula_mod, data = datos, Hess = TRUE)
  
  cat("\n--- ACCIÓN:", var_accion, "---\n")
  print(summary(m))
  return(m)
}


todas_acciones <- c("action_campsite", "action_hut", "action_formalbaths", 
                    "action_formaltrail", "action_trashrecolection", "action_signal", 
                    "action_cuposvisit", "action_nofire", "action_comunication", 
                    "action_information")


resultados_todos <- lapply(todas_acciones, function(x) analizar_aceptabilidad_158(df_encuesta_158, x))




#saco educacion 


library(MASS)
library(dplyr)


analizar_simple_158 <- function(df, var_accion) {
  datos <- df %>% 
    filter(!is.na(!!sym(var_accion))) %>%
    mutate(y = as.factor(!!sym(var_accion)))

  m <- polr(y ~ procedencia + tipo_visitante + edad_encuestado, data = datos, Hess = TRUE)
  

  ctable <- coef(summary(m))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  tabla <- cbind(ctable, "p value" = round(p, 5))
  
  cat("\n--- ACCIÓN:", var_accion, "---\n")
  print(tabla)
  return(m)
}


resultados_limpios <- lapply(todas_acciones, function(x) analizar_simple_158(df_encuesta_158, x))







for(acc in todas_acciones) {
  cat("\n\n********************************************************")
  cat("\nANÁLISIS PARA LA ACCIÓN:", acc)
  cat("\n********************************************************\n")
  
  # Filtramos NAs para esta acción específica
  datos <- df_encuesta_158 %>% 
    filter(!is.na(!!sym(acc))) %>%
    mutate(y = as.factor(!!sym(acc)))
  
  # --- MODELO 1: Procedencia + Tipo de Visitante ---
  cat("\n--- MODELO 1: procedencia + tipo_visitante ---\n")
  try({
    m1 <- polr(y ~ procedencia + tipo_visitante, data = datos, Hess = TRUE)
    ctable1 <- coef(summary(m1))
    p1 <- pnorm(abs(ctable1[, "t value"]), lower.tail = FALSE) * 2
    print(cbind(ctable1, "p value" = round(p1, 5)))
    cat("AIC Modelo 1:", AIC(m1), "\n")
  })
  
  # --- MODELO 2: Procedencia + Edad ---
  cat("\n--- MODELO 2: procedencia + edad_encuestado ---\n")
  try({
    m2 <- polr(y ~ procedencia + edad_encuestado, data = datos, Hess = TRUE)
    ctable2 <- coef(summary(m2))
    p2 <- pnorm(abs(ctable2[, "t value"]), lower.tail = FALSE) * 2
    print(cbind(ctable2, "p value" = round(p2, 5)))
    cat("AIC Modelo 2:", AIC(m2), "\n")
  })
  
  # --- MODELO 3: Solo Procedencia ---
  cat("\n--- MODELO 3: solo procedencia ---\n")
  try({
    m3 <- polr(y ~ procedencia, data = datos, Hess = TRUE)
    ctable3 <- coef(summary(m3))
    p3 <- pnorm(abs(ctable3[, "t value"]), lower.tail = FALSE) * 2
    print(cbind(ctable3, "p value" = round(p3, 5)))
    cat("AIC Modelo 3:", AIC(m3), "\n")
  })
}





todas_acciones <- c("action_campsite", "action_hut", "action_formalbaths", 
                    "action_formaltrail", "action_trashrecolection", "action_signal", 
                    "action_cuposvisit", "action_nofire", "action_comunication", 
                    "action_information")

# Dataframe para guardar los AIC
tabla_comparativa_aic <- data.frame()

for(acc in todas_acciones) {
  
  datos_limpios <- df_encuesta_158 %>% 
    filter(!is.na(!!sym(acc)), 
           !is.na(procedencia), 
           !is.na(tipo_visitante), 
           !is.na(edad_encuestado), 
           !is.na(educacion_encuestado)) %>%
    mutate(y = as.factor(!!sym(acc)))
  
  try({
    # M1: Procedencia + Tipo Visitante
    m1 <- polr(y ~ procedencia + tipo_visitante, data = datos_limpios, Hess = TRUE)
    # M2: Procedencia + Edad
    m2 <- polr(y ~ procedencia + edad_encuestado, data = datos_limpios, Hess = TRUE)
    # M3: Solo Procedencia
    m3 <- polr(y ~ procedencia, data = datos_limpios, Hess = TRUE)
    # M4: Procedencia + Tipo + Edad (Sin Educación)
    m4 <- polr(y ~ procedencia + tipo_visitante + edad_encuestado, data = datos_limpios, Hess = TRUE)
    # M5: Global (Con Educación)
    m5 <- polr(y ~ procedencia + tipo_visitante + edad_encuestado + educacion_encuestado, data = datos_limpios, Hess = TRUE)
    
    # 3. Guardamos los AIC en la tabla
    nueva_fila <- data.frame(
      Accion = acc,
      n_obs = nrow(datos_limpios),
      AIC_M1_ProcTipo = AIC(m1),
      AIC_M2_ProcEdad = AIC(m2),
      AIC_M3_SoloProc = AIC(m3),
      AIC_M4_SinEdu = AIC(m4),
      AIC_M5_Global = AIC(m5)
    )
    tabla_comparativa_aic <- rbind(tabla_comparativa_aic, nueva_fila)
  }, silent = TRUE)
}


print(tabla_comparativa_aic)


#VALIDACION 


library(MASS)
library(DescTools)
library(brant)
library(car)
library(dplyr)
library(performance)


formulas_lista <- list(
  M1_ProcTipo = "y ~ procedencia + tipo_visitante",
  M2_ProcEdad = "y ~ procedencia + edad_encuestado",
  M3_SoloProc = "y ~ procedencia",
  M4_SinEdu   = "y ~ procedencia + tipo_visitante + edad_encuestado",
  M5_Global   = "y ~ procedencia + tipo_visitante + edad_encuestado + educacion_encuestado"
)


for(acc in todas_acciones) {
  cat("\n\n################################################################")
  cat("\n  ANÁLISIS ESTRATÉGICO PARA:", toupper(acc))
  cat("\n################################################################\n")
  
  # Filtro de consistencia: mismas filas para los 5 modelos de esta acción
  df_clean <- df_encuesta_158 %>%
    filter(!is.na(!!sym(acc)), !is.na(procedencia), !is.na(tipo_visitante), 
           !is.na(edad_encuestado), !is.na(educacion_encuestado)) %>%
    mutate(y = as.factor(!!sym(acc)))
  
  # Dataframe temporal para comparar esta acción
  comp_aic <- data.frame()
  
  # 3. Bucle interno por los 5 modelos
  for(f_name in names(formulas_lista)) {
    cat("\n---------------------------------------------------")
    cat("\nMODELO:", f_name)
    cat("\n---------------------------------------------------\n")
    
    try({
      # A. Ajuste
      mod <- polr(as.formula(formulas_lista[[f_name]]), data = df_clean, Hess = TRUE)
      
      # B. Summary con P-values
      ctable <- coef(summary(mod))
      p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
      cat("\n[COEFICIENTES Y P-VALORES]\n")
      print(round(cbind(ctable, "p value" = p), 5))
      
      # C. Validación: R2 y AIC
      r2 <- PseudoR2(mod, which = "McFadden")
      aic_val <- AIC(mod)
      cat("\nAIC:", round(aic_val, 2), "| Pseudo R2 (McFadden):", round(r2, 4), "\n")
      
      # D. Validación: Test de Brant (proporcionalidad)
      cat("\n[TEST DE BRANT]\n")
      try(print(brant(mod)))
      
      # E. Validación: Precisión (Accuracy)
      pred <- predict(mod)
      acc_val <- sum(diag(table(df_clean$y, pred))) / nrow(df_clean)
      cat("\nPrecisión del modelo (Accuracy):", round(acc_val * 100, 2), "%\n")
      
      # Guardar para tabla comparativa final de la acción
      comp_aic <- rbind(comp_aic, data.frame(Modelo = f_name, AIC = aic_val, R2 = r2, Acc = acc_val))
      
    }, silent = FALSE)
  }
  

  cat("\n--- RANKING DE MODELOS PARA", acc, "(Ordenados por AIC) ---\n")
  print(comp_aic[order(comp_aic$AIC), ])
}






#luego de seleccionar que modelos estan validados y tienen mejor AIC imprimo los summary




df_encuesta_158 <- df_encuesta_158 %>%
  mutate(
    procedencia = as.factor(procedencia),
    tipo_visitante = as.factor(tipo_visitante)
  )


print_modelo_final <- function(formula, nombre_accion, data) {
  cat("\n========================================================\n")
  cat("MODELO FINAL PARA:", toupper(nombre_accion), "\n")
  cat("Fórmula:", deparse(formula), "\n")
  cat("========================================================\n")
  
  m <- polr(formula, data = data, Hess = TRUE)
  
  ctable <- coef(summary(m))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  tabla_final <- cbind(ctable, "p value" = round(p, 5))
  
  print(tabla_final)
  cat("\nAIC del modelo:", round(AIC(m), 2), "\n")
  return(m)
}

# --- GRUPO 1: MODELOS ADITIVOS (Procedencia + Tipo) ---
# Campsite, Bath, Trail, Signal

m_campsite <- print_modelo_final(as.factor(action_campsite) ~ procedencia + tipo_visitante, 
                                 "action_campsite", df_encuesta_158)

m_baths <- print_modelo_final(as.factor(action_formalbaths) ~ procedencia + tipo_visitante, 
                              "action_formalbaths", df_encuesta_158)

m_trail <- print_modelo_final(as.factor(action_formaltrail) ~ procedencia + tipo_visitante, 
                              "action_formaltrail", df_encuesta_158)

m_signal <- print_modelo_final(as.factor(action_signal) ~ procedencia + tipo_visitante, 
                               "action_signal", df_encuesta_158)


# --- GRUPO 2: MODELOS SIMPLES (Solo Procedencia) ---
# Hut, Cupos, No Fire, Trash, Comunication, Information

m_hut <- print_modelo_final(as.factor(action_hut) ~ procedencia, 
                            "action_hut", df_encuesta_158)

m_cupos <- print_modelo_final(as.factor(action_cuposvisit) ~ procedencia, 
                              "action_cuposvisit", df_encuesta_158)

m_nofire <- print_modelo_final(as.factor(action_nofire) ~ procedencia, 
                               "action_nofire", df_encuesta_158)

m_trash <- print_modelo_final(as.factor(action_trashrecolection) ~ procedencia, 
                              "action_trashrecolection", df_encuesta_158)

m_comunication <- print_modelo_final(as.factor(action_comunication) ~ procedencia, 
                                     "action_comunication", df_encuesta_158)

m_information <- print_modelo_final(as.factor(action_information) ~ procedencia, 
                                    "action_information", df_encuesta_158)


#-------------------------
#graficar la salida de los modelos 
#----------------------------------


library(ggeffects)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Lista de acciones y sus p-valores 
info_significancia <- data.frame(
  accion = todas_acciones,
  p_val = c(0.00597, 0.85468, 0.91516, 0.19700, 0.10064, 0.00323, 0.39284, 0.04024, 0.05982, 0.07149)
) %>%
  mutate(label_sig = case_when(
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    p_val < 0.1  ~ ".",
    TRUE ~ ""
  ))


df_pred_consolidado <- data.frame()

for(i in 1:nrow(info_significancia)) {
  acc <- info_significancia$accion[i]
  asterisco <- info_significancia$label_sig[i]
  
  m <- resultados_limpios[[i]] 
  
  pred <- predict_response(m, terms = "procedencia")
  df_p <- as.data.frame(pred) %>%
    mutate(
      accion_nombre = case_when(
        acc == "action_campsite" ~ "Ordering Campsites",
        acc == "action_hut" ~ "New Huts",
        acc == "action_formalbaths" ~ "Install Toilets",
        acc == "action_formaltrail" ~ "Order Trails",
        acc == "action_trashrecolection" ~ "Trash Collection",
        acc == "action_signal" ~ "Signage",
        acc == "action_cuposvisit" ~ "Visitor Quotas",
        acc == "action_nofire" ~ "Prohibit Fire",
        acc == "action_comunication" ~ "Emergency Comms",
        acc == "action_information" ~ "Prior Info"
      ),
      # Añadimos el asterisco al nombre para el gráfico
      facet_label = paste0(accion_nombre, " ", asterisco)
    )
  
  df_pred_consolidado <- rbind(df_pred_consolidado, df_p)
}

df_pred_consolidado <- df_pred_consolidado %>%
  mutate(
    response.level = factor(response.level, levels = 1:5, 
                            labels = c("Strongly Agree", "Agree", "Neutral", "Disagree", "Strongly Disagree")),
    x = factor(x, levels = c("local", "no-local"), labels = c("Local", "Visitor"))
  )

ggplot(df_pred_consolidado, aes(x = x, y = predicted, fill = response.level)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~facet_label, ncol = 5) + 
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Acceptability of Management Actions by Origin",
    subtitle = "Significance: ** p<0.01, * p<0.05, . p<0.1 (Based on Ordinal Regression)",
    x = "Origin of Visitor",
    y = "Predicted Probability (%)",
    fill = "Response Level"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 0, face = "bold"),
    panel.spacing = unit(1, "lines")
  )




orden_niveles <- c(
  "Ordering Campsites **", 
  "Signage **", 
  "Prohibit Fire *", 
  "Emergency Comms .", 
  "Prior Info .",
  "Install Toilets ", 
  "New Huts ", 
  "Visitor Quotas ", 
  "Trash Collection ", 
  "Order Trails "
)


df_pred_consolidado <- df_pred_consolidado %>%
  mutate(facet_label = factor(facet_label, levels = orden_niveles))


ggplot(df_pred_consolidado, aes(x = x, y = predicted, fill = response.level)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~facet_label, ncol = 5) + 
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = " ",
    subtitle = "Significance: ** p<0.01, * p<0.05, . p<0.1",
    x = " ",
    y = "Predicted Probability (%)",
    fill = "Response Level"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )




#Determinación de la intensidad de uso 


campamentos_por_intensidad <- df_macro_final_clean %>%
  select(nombre, intensidad) %>% 
  arrange(intensidad, nombre)

print(campamentos_por_intensidad, n = 30)


resumen_intensidad <- df_macro_final_clean %>%
  group_by(intensidad) %>%
  summarise(
    cantidad_sitios = n(),
    campamentos = paste(nombre, collapse = ", "),
    .groups = "drop"
  )

print(resumen_intensidad)


visitantes_por_campamento <- datos_acampe_largo %>%
  filter(!is.na(lugar_limpio), lugar_limpio != "Otro") %>%
  mutate(tamañogrupo = as.numeric(tamañogrupo)) %>%
  group_by(lugar_limpio) %>%
  summarise(
    n_grupos_encuestados = n(),
    total_visitantes     = sum(tamañogrupo, na.rm = TRUE),
    promedio_por_grupo   = round(mean(tamañogrupo, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(total_visitantes))

print(visitantes_por_campamento, n = 20)





#graficar validacion 


usos_resumidos <- df_puntos_final %>%
  group_by(id_camp) %>%
  summarise(
    Fogones = sum(name_2 == "firepits", na.rm = TRUE),
    Baños   = sum(name_2 == "bath", na.rm = TRUE),
    Basura  = sum(name_2 == "littertrash", na.rm = TRUE),
    .groups = "drop"
  )


datos_grafico_validacion <- df_macro_final_clean %>%
  select(id_camp, intensidad) %>% 
  left_join(usos_resumidos, by = "id_camp") %>%
  pivot_longer(cols = c(Fogones, Baños, Basura), 
               names_to = "Tipo_Impacto", 
               values_to = "Cantidad") %>%
  mutate(intensidad = factor(intensidad, 
                             levels = c("bajo", "medio", "alto")))

paleta_impacto <- RColorBrewer::brewer.pal(4, "YlOrRd")

ggplot(datos_grafico_validacion, aes(x = intensidad, y = Cantidad, fill = intensidad)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  facet_wrap(~Tipo_Impacto, scales = "free_y") +
  scale_fill_manual(values = paleta_impacto, name = "Intensidad") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Intensidad de Uso Categorizada",
    y = "Cantidad de elementos registrados (SIG)"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


kruskal.test(Cantidad ~ intensidad, data = filter(datos_grafico_validacion, Tipo_Impacto == "Fogones"))

kruskal.test(Cantidad ~ intensidad, data = filter(datos_grafico_validacion, Tipo_Impacto == "Baños"))

kruskal.test(Cantidad ~ intensidad, data = filter(datos_grafico_validacion, Tipo_Impacto == "Basura"))

