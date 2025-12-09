
get_percent_high_leverage = function(df, model) {
    n = nrow(df)
    p = length(coef(model))
    h_mean = p / n
    # (h_mean = mean(hatvalues(mod_int_pwo))) # same
    # high leverage
    nrow(df[hatvalues(model) > 2*h_mean,])
    nrow(df[hatvalues(model) > 2*h_mean,]) / n
}


get_percent_outliers = function(df, model) {
    n = nrow(df)
    sum(abs(rstandard(model)) > 2) / n
}

get_percent_high_cooks_distance = function(df, model) {
    n = nrow(df)
    # nrow(df[cooks.distance(model) > 4/n,])
    nrow(df[cooks.distance(model) > 4/n,]) / n
}


get_glm_p_scale_rms_error = function(model, df) {
    y_hat <- predict(glm_mod_pwo, type = "response")
    glm_err <- df$det_freq - y_hat
    sqrt(mean(glm_err^2))
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



plot_poly_mod = function(mod, predictor, df = ebird_features_df_pwo) {
    # 1. Range of precipitation over which to plot the curve
    pred_seq <- seq(min(df[predictor]), max(df[predictor]), length.out = 300)
    # 2. Create a new data frame for prediction
    #    Fill all predictors with their mean values EXCEPT ppt_2024_mm
    newdat <- df[ rep(1, 300), ]             # duplicate one row
    newdat[predictor] <- pred_seq            # replace precip with sequence
    # Replace all other predictors with their means
    other_predictors <- names(df)[ names(df) != predictor ]
    for (p in other_predictors) {
        if (is.numeric(df[[p]])) newdat[[p]] <- mean(df[[p]], na.rm = TRUE)
    }
    # 3. Get predicted log(det_freq + 0.1)
    pred_curve <- predict(mod, newdata = newdat)
    # 4. Plot the scatterplot
    plot(log(df$det_freq + 0.1) ~ df[[predictor]],
         pch = 1,
         col = "black",
         main = paste("log(det_freq) vs.", predictor)
    )
    # 5. Add the polynomial prediction curve
    lines(pred_seq, pred_curve, col = "red", lwd = 3)

    plot(df$det_freq ~ df[[predictor]],
         pch = 1,
         col = "black",
         main = paste("det_freq vs.", predictor)
    )
    lines(pred_seq, exp(pred_curve) - 0.1, col = "blue", lwd = 3)
    # all from ChatGPT
}

plot_poly_mods = function(mod1, mod2, mod3, predictor, df = ebird_features_df_pwo) {
    # 1. Range of precipitation over which to plot the curve
    pred_seq <- seq(min(df[predictor]), max(df[predictor]), length.out = 300)
    # 2. Create a new data frame for prediction
    #    Fill all predictors with their mean values EXCEPT ppt_2024_mm
    newdat <- df[ rep(1, 300), ]             # duplicate one row
    newdat[predictor] <- pred_seq            # replace precip with sequence
    # Replace all other predictors with their means
    other_predictors <- names(df)[ names(df) != predictor ]
    for (p in other_predictors) {
        if (is.numeric(df[[p]])) newdat[[p]] <- mean(df[[p]], na.rm = TRUE)
    }
    # 3. Get predicted log(det_freq + 0.1)
    pred_curve_1 <- predict(mod1, newdata = newdat)
    pred_curve_2 <- predict(mod2, newdata = newdat)
    pred_curve_3 <- predict(mod3, newdata = newdat)

    # 4. Plot the scatterplot
    plot(log(df$det_freq + 0.1) ~ df[[predictor]],
         pch = 1,
         col = "black",
         main = paste("log(det_freq) vs.", predictor)
    )
    # 5. Add the polynomial prediction curve
    lines(pred_seq, pred_curve_1, col = "red", lwd = 3)
    lines(pred_seq, pred_curve_2, col = "blue", lwd = 3)
    lines(pred_seq, pred_curve_3, col = "green", lwd = 3)

    plot(df$det_freq ~ df[[predictor]],
         pch = 1,
         col = "black",
         main = paste("det_freq vs.", predictor)
    )
    lines(pred_seq, exp(pred_curve_1) - 0.1, col = "red", lwd = 3)
    lines(pred_seq, exp(pred_curve_2) - 0.1, col = "blue", lwd = 3)
    lines(pred_seq, exp(pred_curve_3) - 0.1, col = "green", lwd = 3)
    # all from ChatGPT
}