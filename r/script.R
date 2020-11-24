setwd("/group/northoff/share/mg/lzc_paper/r/")
library(feather)
library(mediation)

df = read_feather("med_model_1.feather")
# WARNING: THIS IS MED-SPLIT FOR EACH SUBJECT. IT IS DIFFERENT FROM A GLOBAL MED-SPLIT
df = within(df, {
  mf_rest_lh = integer(mf_rest_lh)
  region = factor(region)
})


fit.totaleffect = lm(lzc_movie ~ lzc_rest, df)
summary(fit.totaleffect)
fit.mediator = lm(mf_rest ~ lzc_rest, df)
summary(fit.mediator)
fit.dv = lm(lzc_movie ~ lzc_rest + mf_rest, df)
summary(fit.dv)


results_movie = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest')
summary(results_movie)

# Causal Mediation Analysis
#
# Quasi-Bayesian Confidence Intervals
#
# Estimate 95% CI Lower 95% CI Upper p-value
# ACME              0.215        0.210         0.22  <2e-16 ***
#   ADE               0.743        0.738         0.75  <2e-16 ***
#   Total Effect      0.958        0.956         0.96  <2e-16 ***
#   Prop. Mediated    0.224        0.219         0.23  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Sample Size Used: 103965
# Simulations: 1000

fit.totaleffect = lm(lzc_ret ~ lzc_rest, df)
summary(fit.totaleffect)
fit.mediator = lm(mf_rest ~ lzc_rest, df)
summary(fit.mediator)
fit.dv = lm(lzc_ret ~ lzc_rest + mf_rest, df)
summary(fit.dv)
results_ret = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest')
summary(results_ret)

# Causal Mediation Analysis
#
# Quasi-Bayesian Confidence Intervals
#
#                Estimate 95% CI Lower 95% CI Upper p-value
# ACME            -0.0943      -0.0993        -0.09  <2e-16 ***
# ADE              0.8369       0.8307         0.84  <2e-16 ***
# Total Effect     0.7426       0.7401         0.75  <2e-16 ***
# Prop. Mediated  -0.1270      -0.1337        -0.12  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Sample Size Used: 103965
#
#
# Simulations: 1000


fit.totaleffect = lm(lzc_movie ~ lzc_rest, df)
summary(fit.totaleffect)
# fit.mediator = lm(mf_rest_lh ~ lzc_rest, df)
fit.mediator = glm(mf_rest_lh ~ lzc_rest, data = df, family = binomial)
summary(fit.mediator)
fit.dv = lm(lzc_movie ~ lzc_rest + mf_rest_lh, df)
summary(fit.dv)
results_movie = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest_lh')
summary(results_movie)

fit.totaleffect = lm(lzc_ret ~ lzc_rest, df)
summary(fit.totaleffect)
# fit.mediator = lm(mf_rest_lh ~ lzc_rest, df)
fit.mediator = glm(mf_rest_lh ~ lzc_rest, data = df, family = binomial)
summary(fit.mediator)
fit.dv = lm(lzc_ret ~ lzc_rest + mf_rest_lh, df)
summary(fit.dv)
results_ret = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest_lh')
summary(results_ret)

library(arrow, warn.conflicts = FALSE)
attach(df)
xz = lzc_rest * mf_rest_lh
summary(lm(lzc_movie ~ lzc_rest + mf_rest_lh + xz))

# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.29910 -0.01525  0.00362  0.02160  0.22764
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.177563   0.001965  90.366  < 2e-16 ***
# lzc_rest     0.783394   0.002263 346.234  < 2e-16 ***
# mf_rest_lh  -0.043609   0.009311  -4.684 2.82e-06 ***
# xz           0.086124   0.008941   9.632  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.04424 on 103961 degrees of freedom
# Multiple R-squared:  0.8584,	Adjusted R-squared:  0.8584
# F-statistic: 2.1e+05 on 3 and 103961 DF,  p-value: < 2.2e-16

summary(lm(lzc_ret ~ lzc_rest + mf_rest_lh + xz))

# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.42620 -0.01552  0.00228  0.02091  0.31702
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.281854   0.002121  132.88   <2e-16 ***
# lzc_rest     0.792895   0.002442  324.64   <2e-16 ***
# mf_rest_lh   0.440683   0.010051   43.84   <2e-16 ***
# xz          -0.429723   0.009652  -44.52   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.04775 on 103961 degrees of freedom
# Multiple R-squared:  0.7562,	Adjusted R-squared:  0.7561
# F-statistic: 1.075e+05 on 3 and 103961 DF,  p-value: < 2.2e-16

df = read_feather("med_model_2.feather")
df = within(df, {
  region = factor(region)
})
attach(df)
xz = lzc_rest * mf_rest_lh

summary(lm(lzc_movie ~ lzc_rest + mf_rest_lh + xz))

# Call:
# lm(formula = lzc_movie ~ lzc_rest + mf_rest_lh + xz)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.296109 -0.015448  0.003055  0.020967  0.227089
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.190597   0.001898 100.441  < 2e-16 ***
# lzc_rest    0.766978   0.002181 351.651  < 2e-16 ***
# mf_rest_lh  0.016091   0.007972   2.018   0.0436 *
# xz          0.034563   0.007684   4.498 6.86e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.04359 on 103961 degrees of freedom
# Multiple R-squared:  0.8625,	Adjusted R-squared:  0.8625
# F-statistic: 2.173e+05 on 3 and 103961 DF,  p-value: < 2.2e-16

summary(lm(lzc_ret ~ lzc_rest + mf_rest_lh + xz))

# Call:
# lm(formula = lzc_ret ~ lzc_rest + mf_rest_lh + xz)
#
# Residuals:
#      Min       1Q   Median       3Q      Max
# -0.42428 -0.01552  0.00237  0.02078  0.31574
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.289202   0.002074  139.41   <2e-16 ***
# lzc_rest     0.783580   0.002384  328.63   <2e-16 ***
# mf_rest_lh   0.436521   0.008715   50.09   <2e-16 ***
# xz          -0.423290   0.008400  -50.39   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.04766 on 103961 degrees of freedom
# Multiple R-squared:  0.7572,	Adjusted R-squared:  0.7571
# F-statistic: 1.08e+05 on 3 and 103961 DF,  p-value: < 2.2e-16

setwd("/group/northoff/share/mg/lzc_paper/r/")
library(feather)
library(mediation)
df = read_feather("med_model_3.feather")

fit.totaleffect = lm(lzc_movie ~ lzc_rest, df)
summary(fit.totaleffect)
fit.mediator = lm(mf_rest ~ lzc_rest, df)
summary(fit.mediator)
fit.dv = lm(lzc_movie ~ lzc_rest + mf_rest, df)
summary(fit.dv)
results_movie = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest')
summary(results_movie)

# Call:
# lm(formula = lzc_movie ~ lzc_rest, data = df)
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.139157  0.002083  0.005258  0.013693  0.066097
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.01711    0.01072  -1.596    0.111
# lzc_rest     1.01112    0.01112  90.924   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.03139 on 715 degrees of freedom
# Multiple R-squared:  0.9204,	Adjusted R-squared:  0.9203
# F-statistic:  8267 on 1 and 715 DF,  p-value: < 2.2e-16
#
# Call:
# lm(formula = mf_rest ~ lzc_rest, data = df)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.054047 -0.026450  0.005188  0.025846  0.070351
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.51718    0.01002  -51.63   <2e-16 ***
# lzc_rest     0.68902    0.01039   66.31   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.02933 on 715 degrees of freedom
# Multiple R-squared:  0.8601,	Adjusted R-squared:  0.8599
# F-statistic:  4397 on 1 and 715 DF,  p-value: < 2.2e-16
#
# Call:
# lm(formula = lzc_movie ~ lzc_rest + mf_rest, data = df)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.131513 -0.003368  0.001586  0.017507  0.064511
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.10461    0.02276   4.597 5.07e-06 ***
# lzc_rest     0.84895    0.02903  29.246  < 2e-16 ***
# mf_rest      0.23535    0.03907   6.024 2.73e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.03064 on 714 degrees of freedom
# Multiple R-squared:  0.9242,	Adjusted R-squared:  0.924
# F-statistic:  4356 on 2 and 714 DF,  p-value: < 2.2e-16
#
# Causal Mediation Analysis
#
# Quasi-Bayesian Confidence Intervals
#
#                Estimate 95% CI Lower 95% CI Upper p-value
# ACME              0.162        0.109         0.22  <2e-16 ***
# ADE               0.849        0.789         0.90  <2e-16 ***
# Total Effect      1.011        0.988         1.03  <2e-16 ***
# Prop. Mediated    0.161        0.108         0.21  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 717
# Simulations: 1000

fit.totaleffect = lm(lzc_ret ~ lzc_rest, df)
summary(fit.totaleffect)
fit.mediator = lm(mf_rest ~ lzc_rest, df)
summary(fit.mediator)
fit.dv = lm(lzc_ret ~ lzc_rest + mf_rest, df)
summary(fit.dv)
results_ret = mediate(fit.mediator, fit.dv, treat='lzc_rest', mediator='mf_rest')
summary(results_ret)

# Call:
# lm(formula = lzc_ret ~ lzc_rest, data = df)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.217107 -0.007801 -0.002179  0.018280  0.052679
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.28931    0.01056   27.39   <2e-16 ***
# lzc_rest     0.77919    0.01096   71.12   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.03092 on 715 degrees of freedom
# Multiple R-squared:  0.8762,	Adjusted R-squared:  0.876
# F-statistic:  5058 on 1 and 715 DF,  p-value: < 2.2e-16
#
# Call:
# lm(formula = mf_rest ~ lzc_rest, data = df)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.054047 -0.026450  0.005188  0.025846  0.070351
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.51718    0.01002  -51.63   <2e-16 ***
# lzc_rest     0.68902    0.01039   66.31   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.02933 on 715 degrees of freedom
# Multiple R-squared:  0.8601,	Adjusted R-squared:  0.8599
# F-statistic:  4397 on 1 and 715 DF,  p-value: < 2.2e-16
#
# Call:
# lm(formula = lzc_ret ~ lzc_rest + mf_rest, data = df)
#
# Residuals:
#       Min        1Q    Median        3Q       Max
# -0.204761 -0.003322  0.001554  0.011524  0.063535
#
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.10934    0.02169   5.040  5.9e-07 ***
# lzc_rest     1.01895    0.02767  36.821  < 2e-16 ***
# mf_rest     -0.34798    0.03725  -9.342  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.02921 on 714 degrees of freedom
# Multiple R-squared:  0.8896,	Adjusted R-squared:  0.8893
# F-statistic:  2878 on 2 and 714 DF,  p-value: < 2.2e-16
#
# Causal Mediation Analysis
#
# Quasi-Bayesian Confidence Intervals
#
#                Estimate 95% CI Lower 95% CI Upper p-value
# ACME             -0.241       -0.292        -0.19  <2e-16 ***
# ADE               1.021        0.969         1.07  <2e-16 ***
# Total Effect      0.780        0.758         0.80  <2e-16 ***
# Prop. Mediated   -0.309       -0.375        -0.24  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Sample Size Used: 717
# Simulations: 1000