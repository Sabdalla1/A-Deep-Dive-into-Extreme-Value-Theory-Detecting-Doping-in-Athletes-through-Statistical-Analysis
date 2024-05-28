install.packages("extRemes")
install.packages("ismev")
install.packages("ggplot2")
install.packages("evmix")
install.packages("pid")
install.packages("ReIns")
library(pid)
library(extRemes)
library(evmix)
library(ismev)
library(ggplot2)
library(evir)
library(ReIns)


female <- read.csv("/Users/seifeldineabdalla/Documents/ProjetS2/sport_femmes.csv")
male <- read.csv("/Users/seifeldineabdalla/Documents/ProjetS2/sport_hommes.csv")
dff <- as.numeric(female$REC_ng_mL)
data <- na.omit(dff)
length(data)


#Etudions le cas des femmes
n <- length(data)
n
tail(sort(data), n = 10) #on regarde les plus grandes valeurs

#On fait un Hillplot pour pouvoir trouver le meilleur k, pour l'estimation de xi
dev.new(width=8, height=6) # ajuster taille de fenêtre graphique 
data <- data[data > 0]
hh <- hillplot(data,orderlim = c(2,500), xlab = 'k', y.alpha = FALSE, try.thresh = NULL, main = 'Hill plot')

abline(v=100, lty = 2)
abline(v=150, lty = 2)
abline(v=200,lty = 2)
abline(v=300,lty = 2)
abline(h=0, lty = 4)


#On fait le fitrange pour trouver la meilleur estimation de xi :

fr <- gpd.fitrange(data, umin= quantile(data,1-800/n), umax = quantile(data,1-15/n)) #peut-être le changer
#On fait le GPD :
gpdfit <- gpd.fit(sort(data, decreasing = TRUE), threshold = quantile(data, 1 - 200/n))
gpd.diag(gpdfit) 
xi_mle <- gpdfit$mle[2]
xi_mle
k = 200

#on regarde la log-likelyhood
gpd.profxi(gpdfit, xlow = -0.05 , xup = 1, conf = 0.95, nint = 100)


#Xi de l'estimateur de Hill
hill_estimator <- Hill(data_clean)
xi_hill <- hill_estimator$gamma[200]
xi_hill 


# Fonction pour calculer les quantiles extrêmes
weissman_quantile <- function(data, gamma_hat, p, k) {
  n <- length(data)
  X_k <- sort(data, decreasing = TRUE)[k]
  q_extreme <- X_k * ((k / (n * p)) ^ gamma_hat)
  return(q_extreme)
}



p = 0.001 #1/1000 Quantille 1-1/1000 (0,1%)
q_extreme <- weissman_quantile(data_clean, xi_hill, p, k)
print(q_extreme)

extreme_values <- data_clean[data_clean > q_extreme]
print(extreme_values)
length(extreme_values)
dev.new(width=8, height=6)

plot(data, main = "Données avec seuil de quantile extrême avec p = 0.001", ylab = "Niveaux d'hormones", xlab = "Index")
abline(h = q_extreme, col = "red", lwd = 2, lty = 2)


### Estimateur de Hill trimmed

trimmed_hill <- function(data, k, k0) {
  n <- length(data)
  data_sorted <- sort(data, decreasing = FALSE)
  log_ratio1 <- log(data_sorted[n - k0] / data_sorted[n - k + 1])
  log_ratio2 <- sum(log(data_sorted[(n - k0):(n - k - 1)] / data_sorted[n - k + 1]))
  
  hill_estimator <- (k0 / (k - k0)) * log_ratio1 + (1 / (k - k0)) * log_ratio2
  return(hill_estimator)
}
# Paramètres
k0 <- 50

# Calcul de l'estimateur de Hill ajusté
#xi_trimmed <- trimmed_hill(data, k, k0)
#xi_trimmed
xii_trimmed <- trimmed_hill(data,k,k0)
xii_trimmed

qextremet <- weissman_quantile(data_clean, xii_trimmed, p, k)
qextremet

extreme_values <- data_clean[data_clean > qextremet]
print(extreme_values)
length(extreme_values)

plot(data, main = "Données avec seuil de quantile extrême pour Hill ajusté avec k0 = 50", ylab = "Niveaux d'hormones", xlab = "Index")
abline(h = qextremet, col = "red", lwd = 2, lty = 2)



