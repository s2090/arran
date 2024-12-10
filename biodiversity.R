library(readxl)
library(tidyverse)
library(stringr)
library(ggplot2)
library(betapart)

# load data and standardise
setwd("C:/Users/mia/Documents/Text/Uni/MSc/Professional Skills/Field Course")
dat <- read_xlsx("Arran.xlsx", sheet="Complete")
dat$individualCount <- as.numeric(dat$individualCount)
dat[is.na(dat$individualCount),]$individualCount <- 1
dat <- dat[dat$individualCount != 0,]
dat$habitatType[dat$habitatType == "Broadleaf"] <- "Woodland"
dat$habitatType[dat$habitatType == "Beech_woodland"] <- "Woodland"
dat$habitatType[dat$habitatType == "Bog_transition"] <- "Bog"
dat <- dat[!is.na(dat$habitatType),]

# pivot data for analysis and transform into presence data
trunc <- dat[c("site", "habitatType", "individualCount", names(dat[25:33]))]
dat.pa <- pivot_longer(trunc, cols = c(kingdom, phylum, class, order, family, genus, species), values_to = "taxon", values_drop_na = TRUE)
names(dat.pa)[names(dat.pa) == "name"] <- "rank"

dat.total <- aggregate(individualCount ~ taxon * rank, data = dat.pa, FUN = "sum")
dat.total$individualCount <- as.numeric(1)
names(dat.total)[names(dat.total) == "individualCount"] <- "presence"

# calculate alpha diversity for both sites at all ranks
biodiversity <- data.frame(rank = "a", alpha = 0, site = "both")
biodiversity <- biodiversity[-1,]
for (i in unique(dat.total$rank)) {
  sums <- sum(dat.total[dat.total$rank == i,]$presence)
  biodiversity <- merge(biodiversity, data.frame(rank = i, alpha = as.numeric(sums), site = "Both"), all = TRUE, sort = FALSE)
}

# calculate alpha biodiversity for sites separately at all ranks
dat.sites <- aggregate(individualCount ~ taxon * rank * site, data = dat.pa, FUN = "sum")
dat.sites$individualCount <- as.numeric(1)
names(dat.sites)[names(dat.sites) == "individualCount"] <- "presence"

biodiversity.sites <- data.frame(rank = "a", alpha = 0, site = "a")
biodiversity.sites <- biodiversity.sites[-1,]
for (j in c("North", "South")) {
  for (i in unique(dat.pa$rank)) {
    sums <- sum(dat.sites[dat.sites$rank == i & dat.sites$site == j,]$presence)
    biodiversity.sites <- merge(biodiversity.sites, data.frame(rank = i, alpha = as.numeric(sums), site = j), all = TRUE, sort = FALSE)
  }
}
biodiversity <- merge(biodiversity, biodiversity.sites, all = TRUE, sort = FALSE)

# calculate alpha biodiversity for each combination of site, habitat and rank
dat.pa <- pivot_wider(dat.pa, id_cols = c(site, habitatType, rank), names_from = c(taxon), values_from = individualCount, values_fn = list(individualCount = ~sum(!is.na(.), na.rm = TRUE)))
dat.pa[4:ncol(dat.pa)][!is.na(dat.pa[4:ncol(dat.pa)])] <- as.numeric(1)

alpha <- data.frame(habitat = "a", rank = "a", south = 0, north = 0)
alpha <- alpha[-1,]
counter <- 0
for (i in dat.pa$habitatType) {
  counter <- counter + 1
  j <- dat.pa$rank[counter]
  if (!is.null(dat.pa[dat.pa$site == "South" & dat.pa$habitatType == i & dat.pa$rank == j,])) {
    s <- rowSums(!is.na(dat.pa[dat.pa$site == "South" & dat.pa$habitatType == i & dat.pa$rank == j, 4:ncol(dat.pa)]))
  } else {
    s <- 0
  }
  if (!is.null(dat.pa[dat.pa$site == "North" & dat.pa$habitatType == i & dat.pa$rank == j,])) {
    n <- rowSums(!is.na(dat.pa[dat.pa$site == "North" & dat.pa$habitatType == i & dat.pa$rank == j, 4:ncol(dat.pa)]))
  } else {
    n <- 0
  }
  if (length(s) == 0) {
    df <- data.frame(habitat = i, rank = j, south = 0, north = n)
  } else {
    if (length(n) == 0) {
      df <- data.frame(habitat = i, rank = j, south = s, north = 0)
    } else {
      df <- data.frame(habitat = i, rank = j, south = s, north = n)
    }
  }
  alpha <- merge(alpha, df, all = TRUE, sort = FALSE)
}

# calculate alpha biodiversity for all habitats at each rank
dat.both <- pivot_longer(trunc, cols = c(kingdom, phylum, class, order, family, genus, species), values_to = "taxon", values_drop_na = TRUE)
dat.both <- pivot_wider(dat.both, id_cols = c(habitatType, name), names_from = c(taxon), values_from = individualCount, values_fn = list(individualCount = ~sum(!is.na(.), na.rm = TRUE)))
names(dat.both)[names(dat.both) == "name"] <- "rank"
dat.both[3:ncol(dat.both)][!is.na(dat.both[3:ncol(dat.both)])] <- 1
tot <- data.frame("a", "a", 0)
names(tot) <- c("habitat", "rank", "both")
tot <- tot[-1,]
for (i in unique(dat.both$habitatType)) {
  for (j in unique(dat.both$rank)) {
    summed <- rowSums(dat.both[dat.both$habitatType == i & dat.both$rank == j, 3:ncol(dat.both)], na.rm = TRUE)
    tot <- merge(tot, data.frame(habitat = i, rank = j, both = summed), all = TRUE, sort = FALSE)
  }
}
alpha <- merge(alpha, tot, all = TRUE, sort = FALSE)
alpha.plot <- pivot_longer(alpha, cols = c(south, north, both), names_to = "site", values_to = "alpha")
alpha.plot <- alpha.plot[!is.na(alpha.plot$alpha) & alpha.plot$alpha != 0,]

# calculate beta diversity and sorensen dissimiliarity between sites per habitat
beta <- data.frame(habitat = "a", rank = "a", simpson = 0, sorensen = 0)
beta <- beta[-1,]
for (i in unique(dat.pa$habitatType)) {
  for (j in unique(dat.pa$rank)) {
    if (!(i %in% c("Bog", "Semi-improved grassland", "Heath")) & !(j %in% c("kingdom", "phylum", "class"))) {
      if (dim(dat.both[dat.both$habitatType == i & dat.both$rank == j,])[1] != 0) {
        dat.beta <- dat.pa[dat.pa$habitatType == i & dat.pa$rank == j, !names(dat.pa) %in% c("site", "habitatType", "rank")]
        dat.beta <- Filter(function(x)!all(is.na(x)), dat.beta)
        dat.beta[!is.na(dat.beta)] <- 1
        dat.beta[is.na(dat.beta)] <- 0
        b <- beta.pair(dat.beta)
        if (length(b$beta.sor) != 0) {
          beta <- merge(beta, data.frame(habitat = i, rank = j, simpson = as.numeric(b$beta.sim), sorensen = as.numeric(b$beta.sor)), all = TRUE, sort = FALSE)
        }
      }
    }
  }
}

# calculate beta diversity and sorensen dissimiliarity between sites
dat.sites.pa <- pivot_wider(dat.sites, id_cols = c(site, rank), names_from = c(taxon), values_from = presence)
dat.sites.pa[3:ncol(dat.sites.pa)][!is.na(dat.sites.pa[3:ncol(dat.sites.pa)])] <- 1
beta.sites <- data.frame(rank = "a", simpson = 0, sorensen = 0)
beta.sites <- beta.sites[-1,]
for (j in unique(dat.sites.pa$rank)) {
  dat.beta <- dat.sites.pa[dat.sites.pa$rank == j, !names(dat.sites.pa) %in% c("site", "rank")]
  dat.beta <- Filter(function(x)!all(is.na(x)), dat.beta)
  dat.beta[!is.na(dat.beta)] <- 1
  dat.beta[is.na(dat.beta)] <- 0
  b <- beta.pair(dat.beta)
  if (length(b$beta.sor) != 0) {
    beta.sites <- merge(beta.sites, data.frame(rank = j, simpson = as.numeric(b$beta.sim), sorensen = as.numeric(b$beta.sor)), all = TRUE, sort = FALSE)
  }
}

ggplot(data = biodiversity, aes(x = factor(str_to_title(rank), levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), y = alpha, fill = site)) +
  geom_col(position = position_dodge2(padding = 0.25), width = 0.7) +
  xlab("Taxonomic rank") + ylab("Alpha diversity") + labs(fill = "Site") +
  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

ggplot(data = alpha.plot[alpha.plot$rank %in% c("order", "family", "genus", "species"),], aes(x = str_to_title(site), y = alpha, fill = str_to_title(habitat))) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.25)) +
  facet_wrap(~factor(str_to_sentence(rank), levels = c("Order", "Family", "Genus", "Species"))) +
  xlab("Site") + ylab("Alpha diversity") + labs(fill = "Habitat") +
  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

ggplot(data = beta[beta$rank == "order",], aes(x = habitat, y = sorensen)) +
  geom_col(position = "dodge", width = 0.75) + ylim(c(0,1)) +
  xlab("Habitat") + ylab("Sørensen dissimilarity") +
  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

ggplot(data = beta.sites, aes(x = factor(str_to_title(rank), levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), y = sorensen)) +
  geom_col(position = "dodge", width = 0.75) + ylim(c(0,1)) +
  xlab("Rank") + ylab("Sørensen dissimilarity") +
  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#beta <- pivot_longer(beta, cols = c(simpson, sorensen), names_to = "method", values_to = "beta diversity")
#beta$`beta diversity` <- as.numeric(beta$`beta diversity`)
#ggplot(data = beta[beta$rank == "order",], aes(x = habitat, y = `beta diversity`, fill = habitat)) +
#  geom_col() +
#  facet_wrap(~method)

#ggplot(data = alpha.plot[alpha.plot$rank %in% c("order", "family", "genus", "species"),], aes(x = str_to_title(site), y = alpha, fill = habitat)) +
#  geom_col(position = position_dodge2(preserve = "single", padding = 0.3)) +
#  facet_wrap(~factor(str_to_sentence(rank), levels = c("Order", "Family", "Genus", "Species"))) +
#  xlab("Habitat") + ylab("Alpha diversity") + labs(fill = "Site") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#ggplot(data = alpha.plot[!(alpha.plot$habitat %in% c("Semi-improved grassland", "Heath", "Stream")) & alpha.plot$rank %in% c("order", "family", "genus", "species"),], aes(x = habitat, y = alpha, fill = str_to_title(site))) +
#  geom_col(position = position_dodge2(preserve = "single", padding = 0.3)) +
#  facet_wrap(~factor(str_to_sentence(rank), levels = c("Order", "Family", "Genus", "Species"))) +
#  xlab("Habitat") + ylab("Alpha diversity") + labs(fill = "Site") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#ggplot(data = beta[beta$rank %in% c("order", "family", "genus", "species"),], aes(x = habitat, y = sorensen, fill = habitat)) +
#  geom_col(position = "dodge") +
#  facet_wrap(~factor(str_to_sentence(rank), levels = c("Order", "Family", "Genus", "Species"))) +
#  xlab("Habitat") + ylab("Sørensen dissimilarity") + labs(fill = "Habitat") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#ggplot(data = alpha.plot[alpha.plot$rank == "order",], aes(x = str_to_title(site), y = alpha, fill = habitat)) +
#  geom_col(position = position_dodge2(preserve = "single", padding = 0.3)) +
#  xlab("Habitat") + ylab("Order diversity") + labs(fill = "Site") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#ggplot(data = alpha.plot[alpha.plot$rank == "order",], aes(x = habitat, y = alpha, fill = str_to_title(site))) +
#  geom_col(position = position_dodge2(preserve = "single", padding = 0.3)) +
#  xlab("Habitat") + ylab("Order diversity") + labs(fill = "Site") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))

#ggplot(data = alpha.plot[alpha.plot$rank == "family",], aes(x = habitat, y = alpha, fill = str_to_title(site))) +
#  geom_col(position = position_dodge2(preserve = "single", padding = 0.3)) +
#  xlab("Habitat") + ylab("Family diversity") + labs(fill = "Site") +
#  theme(panel.background = element_rect(fill = "#FAFAFA"), panel.grid = element_line(color = "grey"))