#---------------------- 
rm(list = ls())
require(ggplot2); require(dplyr)

IF = function(u, Nhat, Yhat, Hhat) {
  return(1/Nhat + (log(u)-1)*u/Yhat - Hhat*u/Yhat^2)
}

IF_w1 = function(u, Nhat, Yhat, Hhat) {
  p=-(exp(-Hhat/Yhat-1)*Yhat)/Nhat
  return(-Yhat/(Nhat* pracma::lambertWn(p)))
}

IF_w2 = function(u, Nhat, Yhat, Hhat) {
  p=-(exp(-Hhat/Yhat-1)*Yhat)/Nhat
  return(-Yhat/(Nhat* pracma::lambertWp(p)))
}

StimaTheil <- function(dom) {
  IF.curve = list()
  subdomains = unique(d[[dom]])
  l = 1
  Theilhat = c(); Vhat = c(); W0 = c(); W1 = c()
  for(subdomain in subdomains){ 
    try({
      id.x = d[[dom]]==subdomain
      d_subdomain = d[id.x,]
      
      y = d_subdomain$red_eq
      
      p_i = 1/d$weight
      p_i = p_i[id.x]
      
      w = p_i^(-1)
      Nhat = sum(w) #
      Yhat = sum(y*w)
      Hhat = sum(y*log(y)*w)
      
      IF_y = IF(u = y, Nhat, Yhat, Hhat);  
      ggplot(data.frame(y, IF_y), aes(x = y, y = IF_y)) + geom_line() + geom_point()
      
      Vhat[l] = jipApprox::HTvar(IF_y, pikl = p_ij[id.x, id.x], sample = T, method = 'HT')
      
      Theilhat[l] = dineq::theil.wtd(x = y, weights = w)
      print(sqrt(Vhat[l])/Theilhat[l]*100)
      
      p=-(exp(-Hhat/Yhat-1)*Yhat)/Nhat
      w0 = W0[l] = -Yhat/(Nhat* pracma::lambertWn(p))
      w1 = W1[l] = -Yhat/(Nhat* pracma::lambertWp(p))
      print(w1)
      
      IF.curve[[l]] = data.frame(Domain = subdomain, y = y, influence = IF_y, w0 = w0, w1 = w1); print(IF.curve[[l]])
      l = l+1
    })
  }
  tuscany = c(dineq::theil.wtd(x = d$red_eq, weights = d$weight),
              jipApprox::HTvar(IF(u = d$red_eq, Nhat = sum(d$weight),
                                  Yhat = sum(d$red_eq*d$weight),
                                  Hhat = sum(d$red_eq*log(d$red_eq)*d$weight)), pikl = p_ij, sample = T, method = 'HT'),
              -sum(d$red_eq*d$weight)/(sum(d$weight)* pracma::lambertWn(p)),
              -sum(d$red_eq*d$weight)/(sum(d$weight)* pracma::lambertWp(p))
  )
  
  return(list(Theilhat=Theilhat, Vhat=Vhat,  IF.curve = IF.curve, W0 = W0, W1 = W1, Tuscany = tuscany))
}

convexity = function(y, ymax, ymin) {
  return( (1/y)*log(ymax/ymin) )
}

#################################################################################
##                        Preamble                                            ###
#################################################################################

d <- readxl::read_excel("IRPET_FINE2.xlsx")
d = d[!is.na(d$prov_num),]
d = d[d$anno == "2021",]
d = as.data.frame(lapply(d, rep, d$ncomp))
y = d$red_eq 
hist(y, breaks = 100)

p_i = 1/d$weight # prob inclusione
w = p_i^(-1); 
# p_ij = outer(p_i, p_i); diag(p_ij) = p_i

Nhat = sum(w) 
Yhat = sum(y*w)
Hhat = sum(y*log(y)*w) 


################################################################################
##                            Tuscany                                         ##
################################################################################

bins = 40
ggplot(d, aes(x = red_eq)) + 
  # geom_density(bw = 270, bounds = c(0, Inf)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', color = 'black', bins = bins) +
  theme_bw() + xlab("")
# ggsave(filename = "tuscany_density.eps", device = "eps", width = 10, height = 10)

out <- StimaTheil("prov_num")

out$Tuscany[1] %>% round(5) # theil
sqrt(out$Tuscany[2]) %>% round(5)# SD
sqrt(out$Tuscany[2])*out$Tuscany[1]*100 %>% round(5) # CV
round(out$Tuscany[1]/log(sum(d$weight)) , 5) # REL
round(sqrt(out$Tuscany[2]/log(sum(d$weight)^2)), 5) # Var Rel Theil
out$Tuscany[3] %>% round(5)
out$Tuscany[4] %>% round(5)

convexity(y = sum(d$red_eq), ymax = max(d$red_eq), ymin = min(d$red_eq))*1e6

################################################################################
##                            Provinces                                       ##
################################################################################

out <- StimaTheil("prov_num")

res = data.frame(Domain = names(table(d$prov_num)), 
                 n = as.vector(table(d$prov_num)), 
                 out$Theilhat, SD = sqrt(out$Vhat), 
                 CV = sqrt(out$Vhat)/out$Theilhat, 
                 relT = out$Theilhat/log(tapply(d$weight, INDEX = d$prov_num, sum)),
                 VrelT = sqrt(out$Vhat/log(tapply(d$weight, INDEX = d$prov_num, sum))^2),
                 W0 = out$W0, W1 = out$W1, convexity = sapply( split(d, f = d$prov_num ), function(x) convexity(y = sum(x$red_eq), ymax = max(x$red_eq), ymin = min(x$red_eq) ) )*1e6, row.names = NULL)

res = res %>% inner_join(d %>% select(prov_name, prov_num) %>% unique(), by = c("Domain" = "prov_num")); res$Domain = res$prov_name; res = res[,-11]

kableExtra::kbl(res[,-c(2,10)] %>% arrange(Domain), 
                digits = 5,
                format = "latex",  
                caption = "Theil index estimates (HBS)", 
                col.names = c("Domain", "$\\hat T$", 
                              "$\\sqrt{\\hat{\\text{Var}}[\\hat T]}$", "CV", 
                              "$\\hat T / \\log(\\hat N)$", "$\\sqrt{\\hat{\\text{Var}}[\\hat T / \\log(\\hat N)]}$",
                              "$u_0$", "$u_1$"),
                escape = F, booktabs = TRUE)

res.map = data.frame(Domain = names(table(d$prov_num)), 
                     n = as.vector(table(d$prov_name)),
                     out$Theilhat, SD = sqrt(out$Vhat), CV = sqrt(out$Vhat)/out$Theilhat*100)

### Figures ####

ggplot(d, aes(x = red_eq)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', color = 'black', bins = bins) +
  facet_wrap(~prov_name) + 
  theme_bw() + xlab("") + 
  theme(strip.text.x = element_text(size = 16))
# ggsave(filename = "provinces_density.eps", device = "eps", width = 10, height = 10)

out <- StimaTheil("prov_name")
require(dplyr)
data.frame(do.call(rbind, out$IF.curve)) %>%
  ggplot(aes(x = y, y = influence)) + geom_line(linewidth = .1) + geom_point(size = .7, alpha = 1) + 
  geom_vline(aes(xintercept = w0), linewidth = .2, linetype = "dashed") + 
  geom_vline(aes(xintercept = w1), linewidth = .2, linetype = "dashed") + 
  facet_wrap(~Domain) + theme_bw() + xlab("") + ylab("Influence Function") + 
  theme(strip.text.x = element_text(size = 16))
# ggsave(filename = "prov_influence.eps", device = "eps", width = 10, height = 10)

data.frame(do.call(rbind, out$IF.curve)) %>% 
  ggplot(aes(x = y, y = influence, color = Domain)) + 
  geom_point() + geom_line() + theme_bw() + 
  xlab("") + ylab("Influence Function") + 
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.7), 
        legend.key.size = unit(0.85, "cm"), legend.text = element_text(size = 16))
# ggsave(filename = "prov_influence_all.eps", device = "eps", width = 10, height = 10)

sf = sf::st_read("shapefile_3005c10eec2dca6c85d2fcbfa93c1def/am_com_multipart.shp")
raccordi = readxl::read_excel("raccordo2.xlsx")

map.prov = sf %>% 
  group_by(CODPROV, SIGLA_PROV) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup() %>%
  inner_join(res.map %>%
              dplyr::select(Domain, out.Theilhat), 
            by = c("CODPROV" = "Domain")) 

require(ggplot2)
pal <- c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c")
map.prov = map.prov %>% mutate(clss=case_when(
  out.Theilhat<0.08~"0.08",
  out.Theilhat<0.09~"0.09",
  out.Theilhat<0.10~"0.10",
  out.Theilhat<0.11~"0.11",
  out.Theilhat<0.15~"0.15"
))

ggplot(data = map.prov %>% na.omit()) + 
  geom_sf(aes(fill = clss), color= "white", linewidth = 0.4, show.legend = TRUE) + 
  geom_sf_label(aes(fill = clss, label = SIGLA_PROV), color = "black", show.legend = FALSE) +
  theme_void() +
  scale_fill_manual(name = "Theil Index", values = pal) +
  guides(fill=guide_legend(nrow = 1, title.position = 'top', label.position = 'bottom')) + 
  theme(legend.position = c(0.5,0.025), 
        legend.key.width = unit(3.5, "cm"), legend.key.height = unit(0.2, "cm"), 
        legend.key.spacing = unit(-0.1, "cm"), 
        legend.text = element_text(hjust = 1, vjust = -1.2), legend.title = element_text(vjust = 1.2)) 

ggsave(filename = "mapProvinces.eps", device = "eps", width = 10, height = 10)


#########################################################################################
##                                     Zone                                       #######
#########################################################################################

out <- StimaTheil("zona_num")

res = data.frame(Domain = names(table(d$Zona)),
                 n = as.vector(table(d$Zona)), 
                 out$Theilhat, SD = sqrt(out$Vhat), 
                 CV = sqrt(out$Vhat)/out$Theilhat, 
                 relT = out$Theilhat/log(tapply(d$weight, INDEX = d$Zona, sum)),
                 VrelT = sqrt(out$Vhat/log(tapply(d$weight, INDEX = d$Zona, sum))^2),
                 W0 = out$W0, W1 = out$W1, 
                 convexity = sapply( split(d, f = d$Zona ), function(x) convexity(y = sum(x$red_eq), ymax = max(x$red_eq), ymin = min(x$red_eq) ) )*1e6, row.names = NULL)
res = res[c(2,3,5,1,6,4),]; rownames(res) = NULL


kableExtra::kbl(res[,-c(2,10)], 
                digits = 5,
                format = "latex",  
                caption = "Theil index estimates (HBS)", 
                col.names = c("Domain", "$\\hat T$", 
                              "$\\sqrt{\\hat{\\text{Var}}[\\hat T]}$", "CV", 
                              "$\\hat T / \\log(\\hat N)$", "$\\sqrt{\\hat{\\text{Var}}[\\hat T / \\log(\\hat N)]}$", "$u_0$", "$u_1$"),
                escape = F, booktabs = TRUE)

kableExtra::kbl(res[,c(1,8:10)], 
                digits = 2,
                format = "latex",  
                caption = "Properties of the Influence Function", 
                col.names = c("Domain",
                              "$u_0$", "$u_1$",
                              "Convexity"),
                escape = F, booktabs = TRUE)


### Figures ####

ggplot(d %>% 
         mutate(
           across(Zona,
                  ~factor(., 
                          levels=c("Cities","Made in Italy","Other Industry", "Seaside Tourism", "Agri-Touristic", "North Apennine")))), aes(x = red_eq)) + 
  geom_histogram(aes(y = ..density..), fill = 'white', color = 'black', bins = bins) +
  facet_wrap(~Zona, labeller = label_wrap_gen(width = 20)) + 
  theme_bw() + xlab("") +  theme(strip.text.x = element_text(size = 16))

ggsave(filename = "zones_density.eps", device = "eps", width = 10, height = 10)

out <- StimaTheil("Zona")
data.frame(do.call(rbind, out$IF.curve)) %>%
  ggplot(aes(x = y, y = influence, color = stringr::str_wrap(Domain, 15))) +  
  geom_point(size = 1.5) + geom_line() +
  theme_bw() + xlab("") + ylab("Influence Function") + 
  scale_color_manual(breaks = c("Cities",
                                "Made in Italy",
                                "Other Industry", 
                                "Seaside Tourism", 
                                "Agri-Touristic", 
                                "North Apennine"), values = c("#00798c","#d1495b", "#edae49", "#66a182","#2e4057","#8d96a3")) +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.75), 
        legend.key.width = unit(0.5, "cm"), 
        legend.key.height = unit(1, "cm"), 
        legend.text = element_text(size = 16))
ggsave(filename = "zones_influence_all.eps", device = "eps", width = 10, height = 10)

data.frame(do.call(rbind, out$IF.curve)) %>%
  mutate( across(Domain,
                 ~factor(., 
                         levels=c("Cities",
                                  "Made in Italy",
                                  "Other Industry", 
                                  "Seaside Tourism", 
                                  "Agri-Touristic", 
                                  "North Apennine")))) %>%
  ggplot(aes(x = y, y = influence)) + geom_line(linewidth=.1) + geom_point(size = .7, alpha = 1) + 
  geom_vline(aes(xintercept = w0), linewidth = .2, linetype = "dashed") + 
  geom_vline(aes(xintercept = w1), linewidth = .2, linetype = "dashed") + 
  facet_wrap(~Domain, labeller = label_wrap_gen(width = 20)) + theme_bw() + xlab("") + ylab("Influence Function") + 
  theme(strip.text.x = element_text(size = 16))
ggsave(filename = "zones_influence.eps", device = "eps", width = 10, height = 10)


map.zona = sf %>% 
  inner_join(raccordi, by = c("CODCOM" = "cod_com6")) %>%
  left_join(res %>% mutate(zona6 = as.character(1:6)) %>%
              dplyr::select(out.Theilhat, Domain, zona6),
            by = c("zona6"))

pal <- c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c")
map.zona = map.zona %>% mutate(clss=case_when(
  out.Theilhat<0.08~"0.08",
  out.Theilhat<0.09~"0.09",
  out.Theilhat<0.10~"0.10",
  out.Theilhat<0.11~"0.11",
  out.Theilhat<0.13~"0.15"
))

ggplot(data = map.zona %>% na.omit()) +
  geom_sf(aes(fill = clss), color = "white", linewidth = 0.4, show.legend = TRUE) +
  geom_sf_label(aes(fill = clss, label = zona6), color = "black", show.legend = FALSE) +
  theme_void() +
  scale_fill_manual(name = "Theil Index", values = pal) +
  guides(fill=guide_legend(nrow = 1, title.position = 'top', label.position = 'bottom')) + 
  theme(legend.position = c(0.5,0.025), 
        legend.key.width = unit(3.5, "cm"), legend.key.height = unit(0.2, "cm"), 
        legend.key.spacing = unit(-0.1, "cm"), 
        legend.text = element_text(hjust = 1, vjust = -1.2), legend.title = element_text(vjust = 1.2)) 
# ggsave(filename = "mapZones.eps", device = "eps", width = 10, height = 10)

d %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T, caption = "Equivalent disposable income in Tuscany, Production areas, and Provinces", label = "descr")

d %>% group_by(Zona) %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T)

d %>% group_by(prov_name) %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T)


# decomposition
ineq.2d::theil.2d(as.data.frame(d), total = "red_eq", weights = "weight")

decomp = ineq.2d::theil.2d(as.data.frame(d), total = "red_eq", feature = "prov_name", weights = "weight")
res_decomp = data.frame(Domain = subdomains,
                        Within = as.numeric(unlist(as.vector(decomp[2:11]))),
                        Between = as.vector(unlist(as.vector(decomp[12:21]))))

round(apply(res_decomp[,2:3],2,sum),4)

kableExtra::kbl(res,
                digits = 4, format = "latex", caption = "Decomposition of the Theil index", booktabs = TRUE, label = "tab:app-DecompTuscany")

","#D2F2FF","#C1EAFF","#A9E0FB","#1A87BD")
map.zona = map.zona %>% mutate(clss=case_when(
  out.Theilhat<0.08~"0.08",
  out.Theilhat<0.09~"0.09",
  out.Theilhat<0.10~"0.10",
  out.Theilhat<0.11~"0.11",
  out.Theilhat<0.13~"0.15"
))

ggplot(data = map.zona %>% na.omit()) +
  geom_sf(aes(fill = clss), color = "black", linewidth = 0.4, show.legend = TRUE) +
  geom_sf_label(aes(fill = clss, label = zona6), color = "black", show.legend = FALSE) +
  theme_void() +
  scale_fill_manual(name = "Theil Index", values = pal) +
  guides(fill=guide_legend(nrow = 1, title.position = 'top', label.position = 'bottom')) + 
  theme(legend.position = c(0.5,0.025), 
        legend.key.width = unit(3.5, "cm"), legend.key.height = unit(0.2, "cm"), 
        legend.key.spacing = unit(-0.1, "cm"), 
        legend.text = element_text(hjust = 1, vjust = -1.2), legend.title = element_text(vjust = 1.2)) 
ggsave(filename = "mapZones.eps", device = "eps", width = 10, height = 10)

d %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T, caption = "Equivalent disposable income in Tuscany, Production areas, and Provinces", label = "descr")

d %>% group_by(Zona) %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T)

d %>% group_by(prov_name) %>%
  summarise(n = n(),
            Min. = min(red_eq, na.rm = T),
            Q1 = quantile(red_eq, probs = .25),
            Median = median(red_eq, na.rm = TRUE),
            Mean = mean(red_eq, na.rm = T),
            Q3 = quantile(red_eq, probs = .75),
            Max. = max(red_eq, na.rm = T),
            Std.Dev = sd(red_eq, na.rm = T)) %>%
  kableExtra::kbl(digits = 2, format = "latex", booktabs = T)

# decomposition
ineq.2d::theil.2d(as.data.frame(d), total = "red_eq", weights = "weight")

decomp = ineq.2d::theil.2d(as.data.frame(d), total = "red_eq", feature = "prov_name", weights = "weight")
res_decomp = data.frame(Domain = subdomains,
                        Within = as.numeric(unlist(as.vector(decomp[2:11]))),
                        Between = as.vector(unlist(as.vector(decomp[12:21]))))

round(apply(res_decomp[,2:3],2,sum),4)

kableExtra::kbl(res,
                digits = 4, format = "latex", caption = "Decomposition of the Theil index", booktabs = TRUE, label = "tab:app-DecompTuscany")

