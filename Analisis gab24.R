## Cargamos librerias necesarias
library(dplyr)
library(ggplot2)
library(DHARMa)
library(lme4)
library(emmeans)
library(multcomp)
library(nlme)


## Cargamos los datos
datos <- read.csv("DatosPG (1).csv")

## Generamos una columna con la combinacion de especie y tratamiento
datos$Trat <- paste(datos$Especie, datos$Tratamiento)


## Medidas resumen
resumen <- datos |>
  group_by(Trat) |>
  summarise_at(vars(Raices_debajo,Raiz_mas_larga,N_hojas_verdes),
               c("media"=mean,
                 "DE" = sd))

resumen

## Graficos
numerohojas <- ggplot(datos, aes(x=Dias_desde_armado, y=N_hojas_verdes,
                  color=Trat)) +
  geom_smooth(se=F) +
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  labs(x = "Días desde el armado", y = "Número de hojas verdes") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

numerohojas

ggsave("Imagenes/curva_hojas.png", numerohojas,
       width = 16,
       height = 10,
       units = "cm",
       bg = 'transparent'
)

raizmaslarga <- ggplot(datos, aes(x=Dias_desde_armado, y=Raiz_mas_larga,
                              color=Trat)) +
  geom_smooth(se=F) +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  labs(x = "Días desde el armado", y = "Raíz más larga (cm)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

raizmaslarga

ggsave("Imagenes/raiz_mas_larga.png", raizmaslarga,
       width = 16,
       height = 10,
       units = "cm",
       bg = 'transparent'
)


datos$ID_biorollo <- C(datos$ID_biorollo,"contr.treatment" )

datos$Dias_desde_armado <- C(datos$Dias_desde_armado, "contr.treatment" )

## Modelo poisson
modelo_poisson <- glmer(N_hojas_verdes ~ 1 +
                          Especie +
                          Tratamiento +
                          Dias_desde_armado +
                          Especie:Tratamiento +
                          Especie:Dias_desde_armado +
                          Tratamiento:Dias_desde_armado +
                          Especie:Tratamiento:Dias_desde_armado+
                          (1|ID_biorollo),
                        family = "poisson",
                        nAGQ=1,
                        data = datos)

summary(modelo_poisson)

car::Anova(modelo_poisson, type = 3)

emm_Especie <- emmeans(modelo_poisson,  ~Especie, type = "response")

cld_Especie <- multcomp::cld(emm_Especie, adjust = "tukey",
                          Letters = c("a","b","c","d","e","f","g","h"))

comp_especies <- ggplot(cld_Especie, aes(x = cld_Especie$Especie, y = cld_Especie$rate)) +
  geom_bar(stat = "identity", fill = "#4abe25") +
  geom_text(aes(label=cld_Especie$.group),
            vjust=-0.50,
            size = 5) +
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0,18)) +
  labs(x = "Especies", y = "Número de hojas verdes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


emm_dias <- emmeans(modelo_poisson,  ~Dias_desde_armado, type = "response")

cld_dias <- multcomp::cld(emm_dias, adjust = "tukey",
                             Letters = c("a","b","c","d","e","f","g","h"))

comp_dias <- ggplot(cld_dias, aes(x = cld_dias$Dias_desde_armado,y = cld_dias$rate)) +
  geom_bar(stat = "identity", fill = "#be254a") +
  geom_text(aes(label=cld_dias$.group),
            vjust= -0.5,
            size = 5) +
  scale_y_continuous(breaks = seq(0, 15, by = 5), limits = c(0,11) ) +
  labs(x = "Días desde el armado", y = "Número de hojas verdes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

library(patchwork)

comparacion_p <- comp_especies + comp_dias


ggsave("Imagenes/compa_poisson.png", comparacion_p,
       width = 16,
       height = 10,
       units = "cm",
       bg = 'transparent'
)

# Control de problemas de dispersión
simulationOutput <- simulateResiduals(
  fittedModel = modelo_poisson,
  plot = T, refit = F, n=250)

testDispersion(simulationOutput)


## Modelo cociente --

# Calculamos el cociente
datos_c <- datos  |>
  # Creamo una nueva columna con el valor de N_hojas_verdes iniciales (dia 0)
  group_by(Especie, Tratamiento, ID_biorollo)  |>
  mutate(Hojas_iniciales = N_hojas_verdes[Dias_desde_armado == 0][1])  |>
  ungroup()  |>
  # Creamos la nueva variable
  mutate(Ratio_hojas = N_hojas_verdes / Hojas_iniciales)

datos_c$ID_biorollo <- as.factor(datos_c$ID_biorollo )

datos_c$Dias_desde_armado <- as.factor(datos_c$Dias_desde_armado)


datos_c$ID_biorollo <- C(datos_c$ID_biorollo,"contr.treatment" )

datos_c$Dias_desde_armado <- C(datos_c$Dias_desde_armado, "contr.treatment")

# Filtramos los datos del dia 0 porque no tiene variabilidad en el cocinte
datos_c <- datos_c |> filter(!Dias_desde_armado==0)

# Ajustamos el modelo.
modelo_c <- gls(Ratio_hojas ~ 1 + Especie*Tratamiento*Dias_desde_armado,
                weights = varComb(varExp(form=~fitted(.)|Dias_desde_armado)),
                correlation = corCAR1(form=~as.integer(
                  as.numeric(Dias_desde_armado))|ID_biorollo),
                method = "REML",
                data = datos_c
)


anova(modelo_c)


emm_triple <- emmeans(modelo_c,  ~Especie + Tratamiento +
                        Dias_desde_armado,type = "response")

cld_dias <- multcomp::cld(emm_triple, adjust = "tukey",
                          Letters = c("a","b","c","d","e","f","g","h"))



plot(cld_dias)

# Control Supuesto modelo
plot(modelo_c)

res <- resid(modelo_c,type="pearson")
resiudos <- res/sd(res,na.rm=T)

qqnorm(resiudos, main = "")
qqline(resiudos, col = "red")

df <- as.data.frame(resiudos)

qqplot <- ggplot(df,aes(sample=resiudos))+
  ggplot2::stat_qq()+
  ggplot2::stat_qq_line(color="#4abe25", linewidth = 1 ) +
  ggplot2::coord_equal(ratio = 1) +
  ggplot2::labs(x="Quantiles Teoricos", y="Quantiles Muestrales") +
  ggplot2::theme_minimal() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


disper <- ggplot(df,aes(y=resiudos, x = modelo_c$fitted))+
  geom_point(color="#4abe25") +
  #ggplot2::coord_equal(ratio = 1) +
  ggplot2::labs(x="Valores predichos", y="Residuos estudentizados") +
  ggplot2::theme_minimal() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

supuestos <- qqplot + disper

ggsave("Imagenes/supuestos.png", supuestos,
       width = 16,
       height = 10,
       units = "cm",
       bg = 'transparent'
)

# Cargamos la base con la comparacion de medias por DGC
com_dgc <- read.table("cociente dgc.txt",
                      sep = "\t",
                      header = T)


comparacion_c <- ggplot(com_dgc, aes(x = as.factor(Dias_desde_armado),
                    y =  Medias,
                    color = Especie,
                    linetype = Tratamiento,
                    group = paste(Tratamiento, Especie)
                    )) +
  geom_line(linewidth = 1) +
  geom_text(aes(label = Grupos),
            vjust= -0.5, colour = "#000", check_overlap = T,
            size = 5) +
  scale_color_manual(values = c("#be254a", "#4abe25", "#254abe"))+
  labs(x = "Días desde el armado", y = "Cociente") +
  geom_point(aes(x = 0.5, y = 1), color = "black", size = 4)  +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


ggsave("Imagenes/compa_cociente.png", comparacion_c,
       width = 20,
       height = 15,
       units = "cm",
       bg = 'transparent'
)




