library('lme4')
library(readr)
require(MASS)
require(car)
require("DescTools")
library(emmeans)
library(multcomp)

require(effsize)

# library(data.table)
# name = "outputs/roi_lzc_bp-C_violin_data.rest_movie_ret.7.csv"
# name = "outputs/roi_med_freq_bp-C_violin_data.rest_movie_ret.7.csv"
data = read_csv(name)
data = data[data$Condition %in% c("REST", "MOVIE", "RET"), ]
data = within(data, {
  Condition = relevel(factor(Condition), ref = "REST")
  Network = relevel(factor(Network), ref = "DMN")
})

# DT = data.table(data)
# DT[,
#    .(W = shapiro.test(Value)$statistic, P.value = shapiro.test(Value)$p.value),
#    by = .(Condition, Network)]

rest_dmn = data[data$Condition == "REST" & data$Network == "DMN", ]$Value
rest_vis = data[data$Condition == "REST" & data$Network == "VIS", ]$Value
rest_aud = data[data$Condition == "REST" & data$Network == "AUD", ]$Value
rest_mot = data[data$Condition == "REST" & data$Network == "MOT", ]$Value
mov_dmn = data[data$Condition == "MOVIE" & data$Network == "DMN", ]$Value
mov_vis = data[data$Condition == "MOVIE" & data$Network == "VIS", ]$Value
mov_aud = data[data$Condition == "MOVIE" & data$Network == "AUD", ]$Value
mov_mot = data[data$Condition == "MOVIE" & data$Network == "MOT", ]$Value
ret_dmn = data[data$Condition == "RET" & data$Network == "DMN", ]$Value
ret_vis = data[data$Condition == "RET" & data$Network == "VIS", ]$Value
ret_aud = data[data$Condition == "RET" & data$Network == "AUD", ]$Value
ret_mot = data[data$Condition == "RET" & data$Network == "MOT", ]$Value

wilcox.test(rest_dmn, mov_dmn, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_dmn, mov_dmn, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_dmn, mov_dmn, paired = TRUE, conf.level = 0.99)
wilcox.test(rest_dmn, ret_dmn, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_dmn, ret_dmn, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_dmn, ret_dmn, paired = TRUE, conf.level = 0.99)

wilcox.test(rest_vis, mov_vis, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_vis, mov_vis, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_vis, mov_vis, paired = TRUE, conf.level = 0.99)
wilcox.test(rest_vis, ret_vis, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_vis, ret_vis, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_vis, ret_vis, paired = TRUE, conf.level = 0.99)

wilcox.test(rest_aud, mov_aud, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_aud, mov_aud, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_aud, mov_aud, paired = TRUE, conf.level = 0.99)
wilcox.test(rest_aud, ret_aud, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_aud, ret_aud, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_aud, ret_aud, paired = TRUE, conf.level = 0.99)

wilcox.test(rest_mot, mov_mot, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_mot, mov_mot, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_mot, mov_mot, paired = TRUE, conf.level = 0.99)
wilcox.test(rest_mot, ret_mot, alternative = "two.sided", conf.int = TRUE, conf.level = 0.99, paired = TRUE)
t.test(rest_mot, ret_mot, alternative = "two.sided", conf.level = 0.99, paired = TRUE)
cohen.d(rest_mot, ret_mot, paired = TRUE, conf.level = 0.99)

require(sjstats)
rest = data[data$Condition == "REST", ]
rest$Condition = NULL
rest.aov = aov(Value ~ Network, data = rest)
summary(rest.aov)
summary(glht(rest.aov, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
anova_stats(rest.aov)
kruskal.test(Value ~ Network, data = rest)

# =====================================================================
name = "Network_med_freq_bp-C2_long.rest_movie_ret.7.csv"
data = read_csv(name)
data = within(data, {
  Condition = relevel(factor(Condition), ref = "REST")
  Network = relevel(factor(Network), ref = "DMN")
  Subject = factor(Subject)
  # Value = FisherZ(Value)
})

qqp(data$Value, "lnorm")
qqp(data$Value, "norm")
qqp(data$Value, "lnorm")
nbinom = fitdistr(data$Value, "Negative Binomial")
qqp(data$Value, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
poisson = fitdistr(data$Value, "Poisson")
qqp(data$Value, "pois", poisson$estimate)
gamma = fitdistr(data$Value, "gamma")
qqp(data$Value, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

rest_two = data[data$Condition == "REST" & data$Network %in% c("DMN", "non-DMN"), ]
rest_two = within(rest_two, { Condition = factor(Condition); Network = factor(Network)})
rest_two$Condition = NULL
rest_two_lmm = lmer(Value ~ Network + (1 | Subject), data = rest_two)
summary(rest_two_lmm)
Anova(rest_two_lmm)
summary(glht(rest_two_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))

rest_more = data[data$Condition == "REST" & data$Network != "non-DMN", ]
rest_more = within(rest_more, { Condition = factor(Condition); Network = factor(Network)})
rest_more$Condition = NULL
rest_more_lmm = lmer(Value ~ Network + (1 | Subject), data = rest_more)
summary(rest_more_lmm)
Anova(rest_more_lmm)
summary(glht(rest_more_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))

movie_two = data[data$Condition %in% c("REST", "MOVIE") & data$Network %in% c("DMN", "non-DMN"), ]
movie_two = within(movie_two, { Condition = factor(Condition); Network = factor(Network)})
movie_two_lmm = lmer(Value ~ Network + (1 | Subject) + (Network|Condition), data = movie_two)
summary(movie_two_lmm)
Anova(movie_two_lmm)
summary(glht(movie_two_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
summary(glht(movie_two_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

ret_two = data[data$Condition %in% c("REST", "RET") & data$Network %in% c("DMN", "non-DMN"), ]
ret_two = within(ret_two, { Condition = factor(Condition); Network = factor(Network)})
ret_two_lmm = lmer(Value ~ Network + (1 | Subject) + (Network|Condition), data = ret_two)
summary(ret_two_lmm)
Anova(ret_two_lmm)
summary(glht(ret_two_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))

------------------------------------------------------------------------------------------------------

# library(emmeans)
# emmeans(model, list(pairwise ~ Network), adjust = "tukey")
# library(lmerTest)
# difflsmeans(model, test.effs = "Network")

movie_dmn = data[data$Condition %in% c("REST", "MOVIE") & data$Network == "DMN", ]
movie_dmn$Network = NULL
movie_dmn = within(movie_dmn, { Condition = factor(Condition) })
movie_dmn_lmm = lmer(Value ~ Condition + (1 | Subject), data = movie_dmn)
summary(movie_dmn_lmm)
Anova(movie_dmn_lmm)
summary(glht(movie_dmn_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

ret_dmn = data[data$Condition %in% c("RET", "MOVIE") & data$Network == "DMN", ]
ret_dmn$Network = NULL
ret_dmn = within(ret_dmn, { Condition = factor(Condition) })
ret_dmn_lmm = lmer(Value ~ Condition + (1 | Subject), data = ret_dmn)
summary(ret_dmn_lmm)
Anova(ret_dmn_lmm)
summary(glht(ret_dmn_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

movie_more = data[data$Condition %in% c("REST", "MOVIE") & data$Network != "non-DMN", ]
movie_more$Condition = relevel(movie_more$Condition, ref = "REST")
movie_more$Network = relevel(movie_more$Network, ref = "DMN")
movie_more_lmm = lmer(Value ~ Network + Condition + Network:Condition + (1 | Subject), data = movie_more)
summary(movie_more_lmm)
Anova(movie_more_lmm)
summary(glht(movie_more_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
summary(glht(movie_more_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

ret_more = data[data$Condition %in% c("REST", "RET") & data$Network != "non-DMN", ]
ret_more = within(ret_more, { Condition = factor(Condition); Network = factor(Network)})
ret_more_lmm = lmer(Value ~ Network + Condition + Network*Condition + (1 | Subject), data = ret_more)
summary(ret_more_lmm)
Anova(ret_more_lmm)
summary(glht(ret_more_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
summary(glht(ret_more_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

movie_other = data[data$Condition %in% c("REST", "MOVIE") & data$Network %in% c("PVIS", "AUD", "MOT", "FPC"), ]
movie_more = within(movie_other, { Condition = factor(Condition); Network = factor(Network)})
movie_other_lmm = lmer(Value ~ Network + Condition + Network*Condition + (1 | Subject), data = movie_other)
summary(movie_other_lmm)
Anova(movie_other_lmm)
summary(glht(movie_other_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
summary(glht(movie_other_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))

ret_other = data[data$Condition %in% c("RET", "MOVIE") & data$Network %in% c("PVIS", "AUD", "MOT", "FPC"), ]
ret_other = within(ret_other, { Condition = factor(Condition); Network = factor(Network)})
ret_other_lmm = lmer(Value ~ Network + Condition + Network*Condition + (1 | Subject), data = ret_other)
summary(ret_other_lmm)
Anova(ret_other_lmm)
summary(glht(ret_other_lmm, linfct = mcp(Network = "Tukey")), test = adjusted("bonferroni"))
summary(glht(ret_other_lmm, linfct = mcp(Condition = "Tukey")), test = adjusted("bonferroni"))