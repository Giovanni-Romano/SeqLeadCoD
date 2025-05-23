# Packages ----
options(warning = 1)
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(reshape2) # melt()
    library(mcclust) # comp.psm()
    library(HDInterval) # HDInterval()
    library(salso)
  }
)

cmndargs = commandArgs(trailingOnly = TRUE)

# Load objects ----
TYPE = cmndargs[1] #"Gnedin"
SEED = cmndargs[2]
NAME = paste0(TYPE, SEED)
OUTFOLDER = "output/tHMM/"
## Model output ----
RDS = readRDS(paste0(OUTFOLDER, "res_", NAME, ".RDS"))
C = RDS$output$C
alpha = RDS$output$alpha[-1, ]
## Data ----
ghe = readRDS("data/rds/dataGHE.RDS") %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
  unite("IDshort", c(Country, Sex), remove = FALSE)
## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = 2000:2021
nyears = length(years)
transitions = paste(years[-nyears], years[-1], sep = "-")
IDs = ghe %>% pull(IDshort) %>% unique %>% sort
nIDs = length(IDs)
ncauses = length(levels(ghe$CauseS))
ntransitions = length(transitions)
niter = (RDS$input$ctr$ctr_mcmc$nchain + RDS$input$ctr$ctr_mcmc$nburnin)
nburnin = RDS$input$ctr$ctr_mcmc$nburnin

# Show/Save flags ----
SHOW = TRUE
SAVE = TRUE
COMPUTEPSM = COMPUTEPPE = FALSE
if (SAVE){
  IMGFOLDER = paste0("img/tHMM/", NAME, "/")
  if (!dir.exists(IMGFOLDER)){
    dir.create(IMGFOLDER, recursive = TRUE)
  }
}

# Palette causes ----
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(ghe$CauseS %>% as.character()), decreasing = T) %>% names

# Traceplots ----
## Traceplots number of clusters ----
ncl_trace = apply(C, c(2, 3), function(x) length(unique(x)))
dimnames(ncl_trace) = list("year" = years, "iteration" = 1:niter)
plt_trcncl = ncl_trace %>% melt(value.name = "ncl") %>% 
  filter(iteration > 100) %>%
  ggplot(aes(x = iteration, y = ncl)) +
  geom_line() +
  geom_vline(xintercept = nburnin, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(101, 5000, 10000, 15000, 20000)) +
  facet_wrap(~ year) +
  labs(x = "Iteration", y = "Number of clusters")

if (SHOW) plt_trcncl

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcncl.pdf"
  ggsave(filename = filename,
         plot = plt_trcncl,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Traceplot alpha ----
# Create a vector pasting years[x] with years[x+1] with a dash in between
dimnames(alpha) = list("transition" = transitions, 
                       "iteration" = 1:niter)
plt_trcalpha = alpha %>% melt(value.name = "alpha") %>% 
  ggplot(aes(x = iteration, y = alpha)) +
  geom_line() +
  geom_vline(xintercept = nburnin, linetype = "dashed", color = "red") +
  facet_wrap(~ transition, scales = "fixed") +
  labs(x = "Iteration", y = "Alpha")

if (SHOW) plt_trcalpha

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcalpha.pdf"
  ggsave(filename = filename,
         plot = plt_trcalpha,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

## Smoothed tracplot alpha ----
colorsalpha = colorRampPalette(c("blue", "green", "yellow", "red"))(ntransitions)
names(colorsalpha) = transitions
# Create the base plot
p = alpha %>% 
  melt(value.name = "alpha") %>% 
  ggplot(aes(x = iteration, y = alpha, color = transition)) + 
  geom_smooth(method = "loess", span = .1, se = FALSE, lwd = 0.7, key_glyph = draw_key_rect) + 
  geom_vline(xintercept = nburnin, linetype = "dashed", color = "black") + 
  scale_color_manual("Transition", values = colorsalpha) + 
  guides(color = guide_legend(override.aes = list(fill = colorsalpha, alpha = 1))) + 
  labs(x = "Iteration", y = "Alpha", title = "Smoothed Traceplot of Alpha")

# Create position for labels
{# Get the data from the smooth lines
  smooth_data <- ggplot_build(p)$data[[1]]
  
  # Extract endpoints for each transition
  endpoints <- smooth_data %>%
    group_by(group) %>%
    slice(which.max(x)) %>%
    ungroup()
  
  # Join with original data to get transition values
  # Assuming transitions correspond to the group numbers
  endpoints$transition <- levels(alpha %>% 
                                   melt(value.name = "alpha") %>% 
                                   filter(iteration %in% seq(0, niter, by = 50)) %>% 
                                   pull(transition))[endpoints$group]
  }

# Add colored labels with black text
plt_trcalphasmooth = p + geom_label(data = endpoints, 
                                    aes(x = x + 100, y = y, label = transition, fill = transition),
                                    cex = 3,            # Font size
                                    color = "black",       # Black text color
                                    fontface = "bold",     # Bold text for better visibility
                                    label.size = 0.5,      # Border size of the rectangle
                                    label.padding = unit(0.2, "lines"), # Padding inside the rectangle
                                    hjust = 0,             # Left alignment
                                    show.legend = FALSE) +
  scale_fill_manual(values = colorsalpha) +  # Use the same color scheme for the fills
  scale_x_continuous(limits = c(0, max(smooth_data$x) * 1.15))  # Extend x-axis for labels

if (SHOW) plt_trcalphasmooth

if (SAVE){
  width = 12; height = 8; zoom = 1
  filename = "trcalphasmooth.pdf"
  ggsave(filename = filename,
         plot = plt_trcalphasmooth,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Summary of alphas ----
summary_alpha = apply(alpha[ , -(1:nburnin)], 1, function(x) {
  HDI = hdi(x, credMass = 0.95) %>% unname()
  c("mean" = mean(x), "sd" = sd(x), HDI_l = HDI[1], HDI_u = HDI[2])
}) %>% t() %>% as.data.frame()

if (SHOW) View(summary_alpha %>% round(2))

if (SAVE){
  totxt <- knitr::kable(summary_alpha %>% round(2), format = "simple")  # Or "markdown"
  writeLines(totxt, 
             con = paste0(IMGFOLDER, "summary_alpha.txt"))
}

# Plot summary with mean as points and HDI (black) and sd (gold) as error bars
plt_summaryalpha = summary_alpha %>% 
  rownames_to_column("transition") %>% 
  mutate(transition = factor(transition, levels = transitions)) %>% 
  ggplot(aes(y = transition, x = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = HDI_l, xmax = HDI_u), width = 0.2, lty = 1) +
  geom_errorbar(aes(xmin = mean-sd, xmax = mean+sd), width = 0.2, col = "gold2", lty = 1) +
  coord_flip() +
  labs(y = "Transition", x = "Alpha", title = "Posterior mean with 95% HDI (black) and sd (gold)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

if (SHOW) plt_summaryalpha
if (SAVE) {
  width = 12; height = 8; zoom = 0.75
  filename = "summary_alpha.pdf"
  ggsave(filename = filename,
         plot = plt_summaryalpha,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
}

# Probability of coclustering and PPE ----
if (COMPUTEPSM){
  PSM = apply(C[ , , -(1:nburnin)], 2, function(x) comp.psm(t(x)), simplify = FALSE)
  PSM = lapply(PSM, function(x){
    dimnames(x) = list("ID1" = IDs, "ID2" = IDs)
    x
  })
  names(PSM) = years
  saveRDS(PSM, file = paste0("output/tHMM/PSM_", NAME, ".RDS"))
} else {
  PSM = readRDS("output/tHMM/PSM_20148.RDS")
}
if (COMPUTEPPE){
  PPE = list()
  for (y in 1:nyears){
    cat(y, "\t")
    PPE[[y]] = salso(C[ , y, -(1:nburnin)] %>% t)
    gc()
  }
  saveRDS(PPE, file = paste0("output/tHMM/salso_", RDS$input$ctr$ctr_mcmc$seed, ".RDS"))
} else {
  PPE = readRDS("output/tHMM/salso_20148.RDS") %>% do.call(cbind, .)
}
dimnames(PPE) = list("ID" = IDs, "year" = years)


PPE.ncl = apply(PPE, 2, max)

plt_trcncl + 
  geom_hline(data = data.frame(year = years, PPE.ncl = PPE.ncl),
             mapping = aes(yintercept = PPE.ncl), col = "gold2", lty = 2)


# PSM %>% melt() %>% rename(year = L1, prob = value) %>% mutate(year = as.integer(year)) %>% 
#   left_join(PPE %>% melt(value.name = "ppe1"), by = c("ID1" = "ID", "year")) %>%
#   left_join(PPE%>% melt(value.name = "ppe2"), by = c("ID2" = "ID", "year")) %>% 
#   filter(year == 2010) %>% 
#   ggplot() +
#   geom_raster(aes(x = ID1, y = ID2, fill = prob)) +
#   facet_grid(ppe2 ~ ppe1, scales = "free", space = "free") + 
#   scale_fill_viridis_c("Prob. CC", begin = 0, end = 1) + 
#   theme(axis.ticks = element_blank())
