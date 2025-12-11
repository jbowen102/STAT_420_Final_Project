
# I actually did most of the work in Jupyter notebooks, not R files, so I transcribed
# what I thought was the most important code here. If you really want to see ALL my
# work, you can clone my GitHub repo from https://github.com/jbowen102/STAT_420_Final_Project
# R_auk_full_dataset.ipynb is where I built the code to extract/filter eBird data.
# ebird_modeling.ipynb is where I built the regression models.

# One dataset you'll have to download because it's too big to include in the Coursera ZIP
# is the "Distance to the nearest coast" file from here: https://oceancolor.gsfc.nasa.gov/resources/docs/distfromcoast/



PROJECT_DIR = getwd()
source("feature_construction.R")
source("model_analysis.R")

ebird_features_df_pwo <- read_feature_df("features--pile_woodp_2024_us-sc.csv")

# eBird data structure
auk_extract = get_auk_extract("pileated_woodpecker_20251130.txt")
ebird_zf_df_filtered_pwo = zerofill_and_clean(auk_extract$ebd_only_df, auk_extract$sed_only_df)
str(ebird_zf_df_filtered_pwo)

cor(select(ebird_features_df_pwo, mean_num_observers, mean_start_time_hr, med_duration, mean_dist_km,
                        dist_to_coast_mi, elevation_ft,
                        tmean_2024_C, ppt_2024_mm, det_freq))

get_percent_high_leverage = function(df, model) {
    n = nrow(df)
    p = length(coef(model))
    h_mean = p / n
    nrow(df[hatvalues(model) > 2*h_mean,]) / n
}

get_percent_outliers = function(df, model) {
    n = nrow(df)
    sum(abs(rstandard(model)) > 2) / n
}

get_percent_high_cooks_distance = function(df, model) {
    n = nrow(df)
    nrow(df[cooks.distance(model) > 4/n,]) / n
}

get_k10_cv_rmse = function(formula, df) {
    train_control_k10 <- caret::trainControl(method = "cv", number = 10)
    mod_add_pwo_cv_fit <- caret::train(
        form = formula,
        data = df,
        weights = n_checklists,
        method = "lm",
        trControl = train_control_k10
     )
     mod_add_pwo_cv_fit$results[1,"RMSE"]
}

get_k10_cv_rmse_logy = function(formula, df) {
    train_control_k10 <- caret::trainControl(method = "cv", number = 10, savePredictions = "final")
    mod_add_pwo_logy_cv_fit <- caret::train(
        form = formula,
        data = df,
        weights = n_checklists,
        method = "lm",
        trControl = train_control_k10
    )

    # RMSE on original detection-frequency scale
    preds <- mod_add_pwo_logy_cv_fit$pred
    preds$pred_original <- exp(preds$pred)-0.1
    preds$obs_original  <- exp(preds$obs)-0.1

    rmse_original <- sqrt(mean((preds$obs_original - preds$pred_original)^2))
    rmse_original
    # ChatGPT helped with this.
}


get_k10_cv_rmse_glm = function(formula, df) {
    train_control_k10 <- caret::trainControl(method = "cv", number = 10, savePredictions = "final")
    glm_mod_cv_fit <- caret::train(
        form = formula,
        data = df,
        method = "glm",
        family = binomial(),
        weights = n_checklists,
        trControl = train_control_k10
    )
    # RMSE on original detection-frequency scale
    preds <- glm_mod_cv_fit$pred
    sqrt(mean((preds$obs - preds$pred)^2)) # This is back in probability space already
}

# Additive model with raw detection frequency as response
mod_add_pwo = lm(det_freq ~ mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                        dist_to_coast_mi + elevation_ft + dominant_eco_group +
                        tmean_2024_C + ppt_2024_mm,
                data = ebird_features_df_pwo, weight = n_checklists)
summary(mod_add_pwo)

set.seed(42)
get_k10_cv_rmse(formula(mod_add_pwo), ebird_features_df_pwo)

qqnorm(resid(mod_add_pwo), col = "red")
qqline(resid(mod_add_pwo))
shapiro.test(resid(mod_add_pwo))

plot(resid(mod_add_pwo) ~ fitted(mod_add_pwo), col = "blue")
lmtest::bptest(mod_add_pwo)

# Using log(y+0.1) transformation. Used offset because there are some 0-detection cells.
mod_add_pwo_logy = lm(log(det_freq + 0.1) ~ (mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                        dist_to_coast_mi + elevation_ft + dominant_eco_group +
                        tmean_2024_C + ppt_2024_mm),
                data = ebird_features_df_pwo, weight = n_checklists)
summary(mod_add_pwo_logy)

qqnorm(resid(mod_add_pwo_logy), col = "red")
qqline(resid(mod_add_pwo_logy))
shapiro.test(resid(mod_add_pwo_logy))

plot(resid(mod_add_pwo_logy) ~ fitted(mod_add_pwo_logy), col = "blue")
lmtest::bptest(mod_add_pwo_logy)

# base binomial model
form_pwo_glm_base = det_freq ~ 
                    (mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                     dist_to_coast_mi + elevation_ft + dominant_eco_group +
                     tmean_2024_C + ppt_2024_mm)
mod_pwo_glm_base = glm(form_pwo_glm_base,
                       weights = n_checklists,
                       data = ebird_features_df_pwo,
                       family = binomial)
summary(mod_pwo_glm_base)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_base, ebird_features_df_pwo)

qqnorm(resid(mod_pwo_glm_base), col = "red")
qqline(resid(mod_pwo_glm_base))
shapiro.test(resid(mod_pwo_glm_base))

plot(resid(mod_pwo_glm_base) ~ fitted(mod_pwo_glm_base), col = "blue")
lmtest::bptest(mod_pwo_glm_base)

car::vif(mod_pwo_glm_base)

car::crPlots(mod_pwo_glm_base)

get_percent_high_leverage(ebird_features_df_pwo, mod_pwo_glm_base)
get_percent_outliers(ebird_features_df_pwo, mod_pwo_glm_base)
get_percent_high_cooks_distance(ebird_features_df_pwo, mod_pwo_glm_base)

# full poly model
form_pwo_glm_full_poly = det_freq ~ (mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                                     I(mean_dist_km^2)+I(mean_dist_km^3) +
                                     dist_to_coast_mi +
                                     I(dist_to_coast_mi^2)+I(dist_to_coast_mi^3) +
                                     elevation_ft +
                                     I(elevation_ft^2) +
                                     dominant_eco_group + tmean_2024_C + ppt_2024_mm)
mod_pwo_glm_full_poly = glm(form_pwo_glm_full_poly,
                            weights = n_checklists,
                            data = ebird_features_df_pwo,
                            family = binomial)
summary(mod_pwo_glm_full_poly)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)

n = nrow(ebird_features_df_pwo)
set.seed(42)
best_bic_mod_pwo_glm_poly = step(mod_pwo_glm_full_poly, direction = "backward", k = log(n), trace = FALSE)
summary(best_bic_mod_pwo_glm_poly)
anova(mod_pwo_glm_base, best_bic_mod_pwo_glm_poly)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)

anova(mod_pwo_glm_base, mod_pwo_glm_full_poly)

n = nrow(ebird_features_df_pwo)
best_bic_mod_pwo_glm_poly = step(mod_pwo_glm_full_poly, direction = "backward", k = log(n), trace = FALSE)
summary(best_bic_mod_pwo_glm_poly)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_base, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly), ebird_features_df_pwo)

anova(mod_pwo_glm_base, best_bic_mod_pwo_glm_poly)
anova(best_bic_mod_pwo_glm_poly, mod_pwo_glm_full_poly)

# full poly and interactive model
form_pwo_glm_full_poly_int = det_freq ~ (mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                                     I(mean_dist_km^2)+I(mean_dist_km^3) +
                                     dist_to_coast_mi +
                                     I(dist_to_coast_mi^2)+I(dist_to_coast_mi^3) +
                                     elevation_ft +
                                     I(elevation_ft^2) +
                                     dominant_eco_group + tmean_2024_C + ppt_2024_mm
                                     )^2
mod_pwo_glm_full_poly_int = glm(form_pwo_glm_full_poly_int,
                            weights = n_checklists,
                            data = ebird_features_df_pwo,
                            family = binomial)
summary(mod_pwo_glm_full_poly_int)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_base, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly), ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly_int, ebird_features_df_pwo)


# poly and curated interactive model
form_pwo_glm_poly_cint = det_freq ~ ((mean_num_observers + mean_start_time_hr + med_duration + mean_dist_km +
                                           I(mean_dist_km^2) + I(mean_dist_km^3) )^2 +
                                          (dist_to_coast_mi + I(dist_to_coast_mi^2)+I(dist_to_coast_mi^3) +
                                           dominant_eco_group + tmean_2024_C + ppt_2024_mm )^2 +
                                          (elevation_ft + I(elevation_ft^2) +
                                           dominant_eco_group + tmean_2024_C + ppt_2024_mm )^2)
mod_pwo_glm_poly_cint = glm(form_pwo_glm_poly_cint,
                            weights = n_checklists,
                            data = ebird_features_df_pwo,
                            family = binomial)
summary(mod_pwo_glm_poly_cint)

n = nrow(ebird_features_df_pwo)
best_bic_mod_pwo_glm_poly_cint = step(mod_pwo_glm_poly_cint, direction = "backward", k = log(n), trace = FALSE)
summary(best_bic_mod_pwo_glm_poly_cint)

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_base, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly), ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_poly_cint, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly_cint), ebird_features_df_pwo)


# poly and curated interactive model 3
# Based on what interaction terms seem the most signficant in first curated-interaction model
form_pwo_glm_poly_cint3 = det_freq ~ (mean_num_observers + mean_start_time_hr + med_duration +
                                      mean_dist_km + I(mean_dist_km^2) + I(mean_dist_km^3) +
                                           mean_num_observers:mean_start_time_hr +
                                           mean_num_observers:mean_dist_km +
                                           mean_start_time_hr:mean_dist_km +
                                           mean_start_time_hr:I(mean_dist_km^2) +
                                           mean_start_time_hr:I(mean_dist_km^3) + 
                                           med_duration:mean_dist_km + 
                                           med_duration:I(mean_dist_km^2) +
                                           med_duration:I(mean_dist_km^3) +
                                           dominant_eco_group + 
                                           dist_to_coast_mi + I(dist_to_coast_mi^2) + I(dist_to_coast_mi^3) + 
                                           elevation_ft + I(elevation_ft^2) +
                                           tmean_2024_C +
                                           ppt_2024_mm + 
                                           dist_to_coast_mi:ppt_2024_mm + 
                                           I(dist_to_coast_mi^2):ppt_2024_mm + 
                                           I(dist_to_coast_mi^3):ppt_2024_mm
                                           )
mod_pwo_glm_poly_cint3 = glm(form_pwo_glm_poly_cint3,
                            weights = n_checklists,
                            data = ebird_features_df_pwo,
                            family = binomial)
summary(mod_pwo_glm_poly_cint3)

n = nrow(ebird_features_df_pwo)
best_bic_mod_pwo_glm_poly_cint3 = step(mod_pwo_glm_poly_cint3, direction = "backward", k = log(n), trace = FALSE)
summary(best_bic_mod_pwo_glm_poly_cint3)

setdiff(names(coef(mod_pwo_glm_poly_cint3)), names(coef(best_bic_mod_pwo_glm_poly_cint3)))

set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_base, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_full_poly, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly), ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(form_pwo_glm_poly_cint3, ebird_features_df_pwo)
set.seed(42)
get_k10_cv_rmse_glm(formula(best_bic_mod_pwo_glm_poly_cint3), ebird_features_df_pwo)

anova(best_bic_mod_pwo_glm_poly_cint3, mod_pwo_glm_poly_cint3)

anova(best_bic_mod_pwo_glm_poly, best_bic_mod_pwo_glm_poly_cint3)

qqnorm(resid(best_bic_mod_pwo_glm_poly_cint3), col = "red")
qqline(resid(best_bic_mod_pwo_glm_poly_cint3))

plot(resid(best_bic_mod_pwo_glm_poly_cint3) ~ fitted(best_bic_mod_pwo_glm_poly_cint3), col = "blue")
