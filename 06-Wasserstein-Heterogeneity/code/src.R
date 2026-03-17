# -------------------------------------------------------------------------
# src.R
# Self-contained functions for the synthetic Gaussian illustration
# -------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

# -------------------------------------------------------------------------
# convenience operators and transforms
# -------------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

psi_identity <- function(t) t
psi_sq <- function(t) t^2

make_psi_bounded_sq <- function(c0) {
  force(c0)
  function(t) t^2 / (t^2 + c0^2)
}

# -------------------------------------------------------------------------
# internal helpers for distance-matrix methods
# -------------------------------------------------------------------------
.check_distance_matrix <- function(D, name = deparse(substitute(D)), tol = 1e-8) {
  D <- as.matrix(D)
  
  if (!is.numeric(D)) {
    stop(name, " must be a numeric matrix.")
  }
  if (nrow(D) != ncol(D)) {
    stop(name, " must be square.")
  }
  if (nrow(D) < 2L) {
    stop(name, " must be at least 2 x 2.")
  }
  if (anyNA(D)) {
    stop(name, " contains missing values.")
  }
  if (max(abs(D - t(D))) > tol) {
    stop(name, " must be symmetric.")
  }
  if (min(D) < -tol) {
    stop(name, " has negative entries.")
  }
  if (max(abs(diag(D))) > tol) {
    warning(name, " does not have zero diagonal; the diagonal will be ignored where required.")
  }
  
  storage.mode(D) <- "double"
  D
}

.apply_psi <- function(x, psi) {
  out <- try(psi(x), silent = TRUE)
  if (!inherits(out, "try-error") && is.numeric(out) && length(out) == length(x)) {
    return(as.numeric(out))
  }
  vapply(x, psi, numeric(1))
}

.kernel_matrix <- function(D, psi) {
  z <- .apply_psi(as.vector(D), psi)
  matrix(z, nrow = nrow(D), ncol = ncol(D), byrow = FALSE)
}

.offdiag_row_means <- function(H) {
  n <- nrow(H)
  (rowSums(H) - diag(H)) / (n - 1)
}

.p_value_from_z <- function(z, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(z)),
    less = stats::pnorm(z),
    greater = stats::pnorm(z, lower.tail = FALSE)
  )
}

# -------------------------------------------------------------------------
# main heterogeneity estimator from a distance matrix
# -------------------------------------------------------------------------
wasserstein_heterogeneity <- function(D,
                                      psi = psi_sq,
                                      alpha = 0.05,
                                      D0 = NULL,
                                      alternative = c("two.sided", "less", "greater"),
                                      obs_names = rownames(D),
                                      return_kernel = FALSE) {
  alternative <- match.arg(alternative)
  D <- .check_distance_matrix(D)
  
  n <- nrow(D)
  H <- .kernel_matrix(D, psi)
  upper <- upper.tri(H, diag = FALSE)
  
  U_n <- mean(H[upper])
  ghat <- .offdiag_row_means(H) - U_n
  sigma2_hat <- 4 * sum(ghat^2) / (n - 1)
  
  if (is.null(obs_names)) {
    obs_names <- seq_len(n)
  }
  
  eccentricity_table <- data.frame(
    obs = obs_names,
    ghat = ghat,
    rank = rank(-ghat, ties.method = "first"),
    stringsAsFactors = FALSE
  )
  eccentricity_table <- eccentricity_table[order(eccentricity_table$rank), ]
  rownames(eccentricity_table) <- NULL
  
  if (sigma2_hat > 0) {
    se <- sqrt(sigma2_hat / n)
    zcrit <- stats::qnorm(1 - alpha / 2)
    ci <- c(lower = U_n - zcrit * se, upper = U_n + zcrit * se)
  } else {
    se <- NA_real_
    ci <- c(lower = NA_real_, upper = NA_real_)
    warning(
      "The estimated asymptotic variance is zero. ",
      "The Gaussian approximation, standard error, and confidence interval may be invalid."
    )
  }
  
  z_stat <- NA_real_
  p_value <- NA_real_
  if (!is.null(D0)) {
    if (sigma2_hat > 0) {
      z_stat <- (U_n - D0) / se
      p_value <- .p_value_from_z(z_stat, alternative = alternative)
    } else {
      warning("Cannot compute a z-test when the estimated variance is zero.")
    }
  }
  
  out <- list(
    n = n,
    estimate = U_n,
    ghat = ghat,
    sigma2_hat = sigma2_hat,
    se = se,
    ci = ci,
    z_stat = z_stat,
    p_value = p_value,
    eccentricity_table = eccentricity_table
  )
  
  if (return_kernel) {
    out$H <- H
  }
  
  class(out) <- "wass_heterogeneity"
  out
}

wasserstein_by_group <- function(D,
                                 group,
                                 psi = psi_sq,
                                 alpha = 0.05,
                                 top_k = 5L) {
  D <- .check_distance_matrix(D)
  
  if (length(group) != nrow(D)) {
    stop("length(group) must equal nrow(D).")
  }
  
  group <- as.factor(group)
  lev <- levels(group)
  all_names <- rownames(D)
  if (is.null(all_names)) {
    all_names <- seq_len(nrow(D))
  }
  
  fits <- setNames(vector("list", length(lev)), lev)
  top_eccentric <- setNames(vector("list", length(lev)), lev)
  
  for (g in lev) {
    idx <- which(group == g)
    if (length(idx) < 2L) {
      stop("Each group must contain at least two observations. Group ", g, " does not.")
    }
    fit_g <- wasserstein_heterogeneity(
      D = D[idx, idx, drop = FALSE],
      psi = psi,
      alpha = alpha,
      obs_names = all_names[idx],
      return_kernel = FALSE
    )
    fits[[g]] <- fit_g
    top_eccentric[[g]] <- utils::head(fit_g$eccentricity_table, top_k)
  }
  
  summary_table <- data.frame(
    group = lev,
    n = vapply(fits, function(x) x$n, numeric(1)),
    estimate = vapply(fits, function(x) x$estimate, numeric(1)),
    sigma2_hat = vapply(fits, function(x) x$sigma2_hat, numeric(1)),
    se = vapply(fits, function(x) x$se, numeric(1)),
    ci_lower = vapply(fits, function(x) x$ci[[1]], numeric(1)),
    ci_upper = vapply(fits, function(x) x$ci[[2]], numeric(1)),
    stringsAsFactors = FALSE
  )
  
  list(
    fits = fits,
    summary = summary_table,
    top_eccentric = top_eccentric
  )
}

wasserstein_two_sample <- function(D,
                                   group,
                                   groupA = NULL,
                                   groupB = NULL,
                                   psi = psi_sq,
                                   alpha = 0.05,
                                   delta0 = 0,
                                   alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  D <- .check_distance_matrix(D)
  
  if (length(group) != nrow(D)) {
    stop("length(group) must equal nrow(D).")
  }
  
  group <- as.factor(group)
  lev <- levels(group)
  if (is.null(groupA) || is.null(groupB)) {
    if (length(lev) != 2L) {
      stop("Specify groupA and groupB when there are more than two groups.")
    }
    groupA <- lev[1]
    groupB <- lev[2]
  }
  
  if (!(groupA %in% lev)) stop("groupA not found in group labels.")
  if (!(groupB %in% lev)) stop("groupB not found in group labels.")
  groupA <- as.character(groupA)
  groupB <- as.character(groupB)
  if (groupA == groupB) stop("groupA and groupB must be different.")
  
  idxA <- which(group == groupA)
  idxB <- which(group == groupB)
  if (length(idxA) < 2L || length(idxB) < 2L) {
    stop("Each comparison group must contain at least two observations.")
  }
  
  fitA <- wasserstein_heterogeneity(D[idxA, idxA, drop = FALSE], psi = psi, alpha = alpha)
  fitB <- wasserstein_heterogeneity(D[idxB, idxB, drop = FALSE], psi = psi, alpha = alpha)
  
  delta_hat <- fitA$estimate - fitB$estimate
  se_diff <- sqrt(fitA$sigma2_hat / fitA$n + fitB$sigma2_hat / fitB$n)
  
  if (se_diff > 0) {
    zcrit <- stats::qnorm(1 - alpha / 2)
    ci <- c(lower = delta_hat - zcrit * se_diff, upper = delta_hat + zcrit * se_diff)
    z_stat <- (delta_hat - delta0) / se_diff
    p_value <- .p_value_from_z(z_stat, alternative = alternative)
  } else {
    se_diff <- NA_real_
    ci <- c(lower = NA_real_, upper = NA_real_)
    z_stat <- NA_real_
    p_value <- NA_real_
  }
  
  list(
    groupA = groupA,
    groupB = groupB,
    nA = fitA$n,
    nB = fitB$n,
    estimateA = fitA$estimate,
    estimateB = fitB$estimate,
    delta_hat = delta_hat,
    se = se_diff,
    ci_lower = ci[[1]],
    ci_upper = ci[[2]],
    z_stat = z_stat,
    p_value = p_value
  )
}

pairwise_wasserstein_two_sample <- function(D,
                                            group,
                                            psi = psi_sq,
                                            alpha = 0.05,
                                            p_adjust = "none") {
  D <- .check_distance_matrix(D)
  group <- as.factor(group)
  lev <- levels(group)
  
  if (length(group) != nrow(D)) {
    stop("length(group) must equal nrow(D).")
  }
  if (length(lev) < 2L) {
    stop("At least two groups are required.")
  }
  
  pairs <- utils::combn(lev, 2, simplify = FALSE)
  out <- lapply(pairs, function(ab) {
    fit <- wasserstein_two_sample(
      D = D,
      group = group,
      groupA = ab[1],
      groupB = ab[2],
      psi = psi,
      alpha = alpha
    )
    data.frame(
      groupA = fit$groupA,
      groupB = fit$groupB,
      nA = fit$nA,
      nB = fit$nB,
      estimateA = fit$estimateA,
      estimateB = fit$estimateB,
      delta_hat = fit$delta_hat,
      se = fit$se,
      ci_lower = fit$ci_lower,
      ci_upper = fit$ci_upper,
      z_stat = fit$z_stat,
      p_value = fit$p_value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$p_value_adjusted <- stats::p.adjust(out$p_value, method = p_adjust)
  rownames(out) <- NULL
  out
}

plugin_stability_bound <- function(approx_error, L) {
  if (!is.numeric(approx_error) || anyNA(approx_error)) {
    stop("approx_error must be a numeric vector of W_2(hat(mu_i), mu_i) values.")
  }
  if (any(approx_error < 0)) {
    stop("approx_error must contain nonnegative values.")
  }
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 0) {
    stop("L must be one nonnegative number.")
  }
  2 * L * mean(approx_error)
}

# -------------------------------------------------------------------------
# synthetic 2D Gaussian measures
# -------------------------------------------------------------------------
gamma_shape_scale_from_mean_sd <- function(mean, sd) {
  if (any(mean <= 0) || any(sd <= 0)) {
    stop("Gamma means and standard deviations must be strictly positive.")
  }
  shape <- (mean / sd)^2
  scale <- sd^2 / mean
  list(shape = shape, scale = scale)
}

rgamma_mean_sd <- function(n, mean, sd) {
  pars <- gamma_shape_scale_from_mean_sd(mean, sd)
  stats::rgamma(n, shape = pars$shape, scale = pars$scale)
}

safe_sd <- function(x) {
  out <- stats::sd(x)
  if (is.na(out) || out <= 1e-8) 1e-8 else out
}

trace_from_variances <- function(v) {
  sum(v)
}

ellipse_coords <- function(mean_vec, sd_vec, level = 2, n = 200L) {
  theta <- seq(0, 2 * pi, length.out = n)
  cbind(
    x = mean_vec[1] + level * sd_vec[1] * cos(theta),
    y = mean_vec[2] + level * sd_vec[2] * sin(theta)
  )
}

default_measure_specs <- function() {
  list(
    A = list(
      components = list(
        list(
          weight = 1.00,
          label = "baseline",
          mean_m = c(-2.5, 0.0),
          sd_m = c(0.15, 0.10),
          mean_s = c(0.55, 0.35),
          sd_s = c(0.05, 0.03)
        )
      )
    ),
    B = list(
      components = list(
        list(
          weight = 1.00,
          label = "diffuse",
          mean_m = c(2.5, 0.0),
          sd_m = c(0.55, 0.35),
          mean_s = c(0.80, 0.50),
          sd_s = c(0.12, 0.08)
        )
      )
    ),
    C = list(
      components = list(
        list(
          weight = 0.85,
          label = "main",
          mean_m = c(0.0, 2.4),
          sd_m = c(0.22, 0.18),
          mean_s = c(0.50, 0.40),
          sd_s = c(0.06, 0.05)
        ),
        list(
          weight = 0.15,
          label = "atypical",
          mean_m = c(1.4, 3.6),
          sd_m = c(0.10, 0.10),
          mean_s = c(1.10, 0.18),
          sd_s = c(0.05, 0.03)
        )
      )
    )
  )
}

true_D_sq_group <- function(group_spec) {
  comps <- group_spec$components
  w <- vapply(comps, function(cc) cc$weight, numeric(1))
  w <- w / sum(w)
  
  means <- do.call(rbind, lapply(comps, function(cc) c(cc$mean_m, cc$mean_s)))
  vars <- do.call(rbind, lapply(comps, function(cc) c(cc$sd_m^2, cc$sd_s^2)))
  
  overall_mean <- colSums(means * w)
  overall_second <- colSums((vars + means^2) * w)
  overall_var <- overall_second - overall_mean^2
  
  2 * trace_from_variances(overall_var)
}

true_D_sq_all_groups <- function(specs) {
  vapply(specs, true_D_sq_group, numeric(1))
}

simulate_component_measures <- function(n, component, group_name) {
  data.frame(
    group = group_name,
    component = component$label %||% "component",
    m1 = stats::rnorm(n, mean = component$mean_m[1], sd = component$sd_m[1]),
    m2 = stats::rnorm(n, mean = component$mean_m[2], sd = component$sd_m[2]),
    s1 = rgamma_mean_sd(n, mean = component$mean_s[1], sd = component$sd_s[1]),
    s2 = rgamma_mean_sd(n, mean = component$mean_s[2], sd = component$sd_s[2]),
    stringsAsFactors = FALSE
  )
}

simulate_group_measures <- function(n, group_spec, group_name) {
  comps <- group_spec$components
  w <- vapply(comps, function(cc) cc$weight, numeric(1))
  w <- w / sum(w)
  comp_idx <- sample(seq_along(comps), size = n, replace = TRUE, prob = w)
  
  out <- vector("list", length(comps))
  keep <- logical(length(comps))
  for (k in seq_along(comps)) {
    nk <- sum(comp_idx == k)
    if (nk > 0L) {
      out[[k]] <- simulate_component_measures(nk, comps[[k]], group_name)
      keep[k] <- TRUE
    }
  }
  out <- do.call(rbind, out[keep])
  out <- out[sample.int(nrow(out)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

w2_matrix_from_parameter_columns <- function(df, mean_cols, sd_cols, ids = NULL) {
  Z <- as.matrix(df[, c(mean_cols, sd_cols), drop = FALSE])
  D <- as.matrix(stats::dist(Z, method = "euclidean"))
  ids <- ids %||% rownames(df) %||% seq_len(nrow(df))
  rownames(D) <- ids
  colnames(D) <- ids
  D
}

w2_matrix_axis_aligned_gaussian <- function(meta) {
  w2_matrix_from_parameter_columns(
    meta,
    mean_cols = c("m1", "m2"),
    sd_cols = c("s1", "s2"),
    ids = meta$id
  )
}

simulate_measure_dataset <- function(seed = 20260307,
                                     n_A = 60L,
                                     n_B = 60L,
                                     n_C = 60L,
                                     specs = default_measure_specs()) {
  set.seed(seed)
  
  group_sizes <- c(A = n_A, B = n_B, C = n_C)
  meta_list <- lapply(names(group_sizes), function(g) {
    simulate_group_measures(group_sizes[[g]], specs[[g]], g)
  })
  names(meta_list) <- names(group_sizes)
  
  meta <- do.call(rbind, meta_list)
  rownames(meta) <- NULL
  obs_in_group <- ave(seq_len(nrow(meta)), meta$group, FUN = seq_along)
  meta$id <- paste0(meta$group, sprintf("%03d", obs_in_group))
  meta <- meta[, c("id", "group", "component", "m1", "m2", "s1", "s2")]
  
  D <- w2_matrix_axis_aligned_gaussian(meta)
  true_D_sq <- true_D_sq_all_groups(specs)
  
  list(
    meta = meta,
    D = D,
    specs = specs,
    true_D_sq = true_D_sq,
    seed = seed,
    group_sizes = group_sizes
  )
}

approximate_measures_by_sampling <- function(meta, inner_n = 50L, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(meta)
  approx <- data.frame(
    id = meta$id,
    group = meta$group,
    component = meta$component,
    m1_hat = numeric(n),
    m2_hat = numeric(n),
    s1_hat = numeric(n),
    s2_hat = numeric(n),
    stringsAsFactors = FALSE
  )
  approx_error <- numeric(n)
  
  for (i in seq_len(n)) {
    Xi <- cbind(
      stats::rnorm(inner_n, mean = meta$m1[i], sd = meta$s1[i]),
      stats::rnorm(inner_n, mean = meta$m2[i], sd = meta$s2[i])
    )
    mhat <- colMeans(Xi)
    shat <- apply(Xi, 2, safe_sd)
    
    approx$m1_hat[i] <- mhat[1]
    approx$m2_hat[i] <- mhat[2]
    approx$s1_hat[i] <- shat[1]
    approx$s2_hat[i] <- shat[2]
    
    approx_error[i] <- sqrt(
      sum((mhat - c(meta$m1[i], meta$m2[i]))^2) +
        sum((shat - c(meta$s1[i], meta$s2[i]))^2)
    )
  }
  
  list(
    approx = approx,
    approx_error = stats::setNames(approx_error, meta$id)
  )
}

w2_matrix_fitted_gaussian <- function(approx) {
  w2_matrix_from_parameter_columns(
    approx,
    mean_cols = c("m1_hat", "m2_hat"),
    sd_cols = c("s1_hat", "s2_hat"),
    ids = approx$id
  )
}

plugin_stability_curve <- function(sim,
                                   inner_grid = c(10L, 25L, 50L, 100L, 200L),
                                   n_rep = 25L,
                                   seed = sim$seed + 1000L) {
  fit_true <- wasserstein_heterogeneity(sim$D, psi = psi_identity)
  
  out <- vector("list", length(inner_grid) * n_rep)
  counter <- 0L
  for (m in inner_grid) {
    for (r in seq_len(n_rep)) {
      counter <- counter + 1L
      ap <- approximate_measures_by_sampling(
        sim$meta,
        inner_n = m,
        seed = seed + 1000L * m + r
      )
      D_hat <- w2_matrix_fitted_gaussian(ap$approx)
      fit_hat <- wasserstein_heterogeneity(D_hat, psi = psi_identity)
      actual_error <- abs(fit_hat$estimate - fit_true$estimate)
      upper_bound <- plugin_stability_bound(ap$approx_error, L = 1)
      out[[counter]] <- data.frame(
        inner_n = m,
        rep = r,
        estimate_true = fit_true$estimate,
        estimate_hat = fit_hat$estimate,
        actual_error = actual_error,
        upper_bound = upper_bound,
        bound_holds = actual_error <= upper_bound + 1e-12,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, out)
  
  summary <- stats::aggregate(
    cbind(actual_error, upper_bound, bound_holds) ~ inner_n,
    data = out,
    FUN = mean
  )
  summary$bound_holds <- as.numeric(summary$bound_holds)
  
  list(raw = out, summary = summary)
}

balanced_degenerate_dataset <- function(m = 30L,
                                        a = 1.5,
                                        sd_vec = c(0.55, 0.35)) {
  meta <- data.frame(
    id = c(paste0("L", sprintf("%03d", seq_len(m))),
           paste0("R", sprintf("%03d", seq_len(m)))),
    group = rep("degenerate", 2L * m),
    component = c(rep("left", m), rep("right", m)),
    m1 = c(rep(-a, m), rep(a, m)),
    m2 = 0,
    s1 = sd_vec[1],
    s2 = sd_vec[2],
    stringsAsFactors = FALSE
  )
  D <- w2_matrix_axis_aligned_gaussian(meta)
  list(meta = meta, D = D)
}

prepare_single_dataset_plot_data <- function(sim, by_sq, transform_table, plugin_curve) {
  meta <- sim$meta
  group_levels <- c("A", "B", "C")
  
  locations_df <- meta
  locations_df$group <- factor(locations_df$group, levels = group_levels)
  locations_df$component <- factor(locations_df$component,
                                   levels = c("baseline", "diffuse", "main", "atypical"))
  locations_df$size_proxy <- 0.6 + 1.2 * ((meta$s1 + meta$s2) / max(meta$s1 + meta$s2))
  
  summary_sq <- by_sq$summary
  summary_sq <- summary_sq[match(group_levels, summary_sq$group), ]
  summary_sq$group <- factor(summary_sq$group, levels = group_levels)
  summary_sq$true_D_sq <- as.numeric(sim$true_D_sq[as.character(summary_sq$group)])
  
  fit_C <- by_sq$fits[["C"]]
  meta_C <- meta[meta$group == "C", , drop = FALSE]
  meta_C$obs_index <- seq_len(nrow(meta_C))
  meta_C$ghat <- as.numeric(fit_C$ghat[meta_C$id])
  meta_C$type <- ifelse(meta_C$component == "atypical", "atypical component", "main component")
  
  ord_abs <- order(abs(meta_C$ghat))
  typical_ids <- meta_C$id[utils::head(ord_abs, 3L)]
  eccentric_ids <- utils::head(by_sq$top_eccentric[["C"]]$obs, 3L)
  show_ids <- c(typical_ids, eccentric_ids)
  show_type <- c(rep("typical", length(typical_ids)), rep("high eccentricity", length(eccentric_ids)))
  
  ellipse_list <- vector("list", length(show_ids))
  label_df <- meta[match(show_ids, meta$id), c("id", "m1", "m2", "s1", "s2", "component")]
  label_df$type <- show_type
  for (i in seq_along(show_ids)) {
    idx <- match(show_ids[i], meta$id)
    E <- ellipse_coords(c(meta$m1[idx], meta$m2[idx]), c(meta$s1[idx], meta$s2[idx]), level = 2)
    ellipse_list[[i]] <- data.frame(
      id = show_ids[i],
      type = show_type[i],
      x = E[, 1],
      y = E[, 2],
      stringsAsFactors = FALSE
    )
  }
  ellipse_df <- do.call(rbind, ellipse_list)
  ellipse_df$type <- factor(ellipse_df$type, levels = c("typical", "high eccentricity"))
  label_df$type <- factor(label_df$type, levels = c("typical", "high eccentricity"))
  
  transform_long <- rbind(
    data.frame(group = transform_table$group,
               transform = "squared",
               estimate = transform_table$squared,
               stringsAsFactors = FALSE),
    data.frame(group = transform_table$group,
               transform = "bounded",
               estimate = transform_table$bounded,
               stringsAsFactors = FALSE)
  )
  transform_long$group <- factor(transform_long$group, levels = group_levels)
  transform_long$transform <- factor(transform_long$transform, levels = c("squared", "bounded"))
  
  plugin_df <- plugin_curve$summary
  plugin_long <- rbind(
    data.frame(inner_n = plugin_df$inner_n,
               series = "mean upper bound",
               value = plugin_df$upper_bound,
               stringsAsFactors = FALSE),
    data.frame(inner_n = plugin_df$inner_n,
               series = "mean actual error",
               value = plugin_df$actual_error,
               stringsAsFactors = FALSE)
  )
  plugin_long$series <- factor(plugin_long$series,
                               levels = c("mean upper bound", "mean actual error"))
  
  list(
    locations = locations_df,
    group_summary = summary_sq,
    groupC_eccentricities = meta_C,
    representative_ellipses = ellipse_df,
    representative_labels = label_df,
    transform_comparison = transform_long,
    plugin_curve = plugin_df,
    plugin_curve_long = plugin_long
  )
}

run_single_dataset_demo <- function(seed = 20260307,
                                    n_A = 60L,
                                    n_B = 60L,
                                    n_C = 60L,
                                    alpha = 0.05,
                                    top_k = 6L,
                                    make_base_plots = FALSE,
                                    save_plot_data = TRUE,
                                    output_dir = NULL) {
  sim <- simulate_measure_dataset(seed = seed, n_A = n_A, n_B = n_B, n_C = n_C)
  
  by_sq <- wasserstein_by_group(
    sim$D,
    group = sim$meta$group,
    psi = psi_sq,
    alpha = alpha,
    top_k = top_k
  )
  
  summary_sq <- by_sq$summary
  summary_sq <- summary_sq[match(c("A", "B", "C"), summary_sq$group), ]
  summary_sq$true_D_sq <- as.numeric(sim$true_D_sq[summary_sq$group])
  summary_sq$error <- summary_sq$estimate - summary_sq$true_D_sq
  summary_sq$abs_error <- abs(summary_sq$error)
  summary_sq$ci_covers_true <- with(summary_sq, ci_lower <= true_D_sq & true_D_sq <= ci_upper)
  
  fits_zero <- lapply(c("A", "B", "C"), function(g) {
    idx <- sim$meta$group == g
    wasserstein_heterogeneity(sim$D[idx, idx, drop = FALSE], psi = psi_sq, D0 = 0)
  })
  names(fits_zero) <- c("A", "B", "C")
  one_sample_table <- data.frame(
    group = names(fits_zero),
    estimate = vapply(fits_zero, function(x) x$estimate, numeric(1)),
    se = vapply(fits_zero, function(x) x$se, numeric(1)),
    z_stat = vapply(fits_zero, function(x) x$z_stat, numeric(1)),
    p_value = vapply(fits_zero, function(x) x$p_value, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  pairwise_sq <- pairwise_wasserstein_two_sample(
    sim$D,
    group = sim$meta$group,
    psi = psi_sq,
    alpha = alpha,
    p_adjust = "holm"
  )
  pairwise_sq$true_delta <- mapply(
    function(a, b) sim$true_D_sq[[a]] - sim$true_D_sq[[b]],
    pairwise_sq$groupA,
    pairwise_sq$groupB
  )
  pairwise_sq$ci_covers_true <- with(pairwise_sq, ci_lower <= true_delta & true_delta <= ci_upper)
  
  upper <- upper.tri(sim$D, diag = FALSE)
  c0 <- as.numeric(stats::median(sim$D[upper]))
  by_bound <- wasserstein_by_group(
    sim$D,
    sim$meta$group,
    psi = make_psi_bounded_sq(c0),
    alpha = alpha,
    top_k = top_k
  )
  
  transform_table <- data.frame(
    group = c("A", "B", "C"),
    squared = by_sq$summary$estimate[match(c("A", "B", "C"), by_sq$summary$group)],
    bounded = by_bound$summary$estimate[match(c("A", "B", "C"), by_bound$summary$group)],
    stringsAsFactors = FALSE
  )
  
  plugin_curve <- plugin_stability_curve(
    sim,
    inner_grid = c(10L, 25L, 50L, 100L, 200L),
    n_rep = 20L,
    seed = seed + 5000L
  )
  
  deg <- balanced_degenerate_dataset(m = 30L, a = 1.5, sd_vec = c(0.55, 0.35))
  fit_deg <- suppressWarnings(
    wasserstein_heterogeneity(deg$D, psi = psi_sq, D0 = 0)
  )
  degenerate_table <- data.frame(
    estimate = fit_deg$estimate,
    sigma2_hat = fit_deg$sigma2_hat,
    se = fit_deg$se,
    stringsAsFactors = FALSE
  )
  
  top_C <- by_sq$top_eccentric[["C"]]
  top_C$component <- sim$meta$component[match(top_C$obs, sim$meta$id)]
  top_C$m1 <- sim$meta$m1[match(top_C$obs, sim$meta$id)]
  top_C$m2 <- sim$meta$m2[match(top_C$obs, sim$meta$id)]
  top_C$s1 <- sim$meta$s1[match(top_C$obs, sim$meta$id)]
  top_C$s2 <- sim$meta$s2[match(top_C$obs, sim$meta$id)]
  prop_atypical_topC <- mean(top_C$component == "atypical")
  
  plot_data <- prepare_single_dataset_plot_data(sim, by_sq, transform_table, plugin_curve)
  
  invisible(list(
    sim = sim,
    by_sq = by_sq,
    summary_sq = summary_sq,
    one_sample_table = one_sample_table,
    pairwise_sq = pairwise_sq,
    transform_table = transform_table,
    top_C = top_C,
    prop_atypical_topC = prop_atypical_topC,
    plugin_curve = plugin_curve,
    degenerate_table = degenerate_table,
    plot_data = plot_data
  ))
}

# -------------------------------------------------------------------------
# ggplot builders
# -------------------------------------------------------------------------
make_2d_measure_ggplots <- function(plot_data) {
  loc <- plot_data$locations
  ell <- plot_data$representative_ellipses
  lab <- plot_data$representative_labels
  summ <- plot_data$group_summary
  ecc <- plot_data$groupC_eccentricities
  trans <- plot_data$transform_comparison
  plugin <- plot_data$plugin_curve_long
  
  group_cols <- c(A = "steelblue4", B = "darkorange3", C = "darkgreen")
  type_cols <- c("typical" = "grey45", "high eccentricity" = "firebrick")
  comp_cols <- c("main component" = "grey35", "atypical component" = "firebrick")
  trans_fills <- c("squared" = "grey35", "bounded" = "grey75")
  plugin_cols <- c("mean upper bound" = "black", "mean actual error" = "grey35")
  plugin_shapes <- c("mean upper bound" = 17, "mean actual error" = 19)
  
  p_locations <- ggplot(loc, aes(x = m1, y = m2, colour = group, size = size_proxy)) +
    geom_point(alpha = 0.95) +
    scale_colour_manual(values = group_cols) +
    scale_size_identity() +
    labs(
      x = expression(m[1]),
      y = expression(m[2]),
      title = "Locations of the generated measures",
      colour = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  p_representative <- ggplot() +
    geom_path(data = ell,
              aes(x = x, y = y, group = id, colour = type),
              linewidth = 0.8) +
    geom_text(data = lab,
              aes(x = m1, y = m2, label = id),
              vjust = -0.4,
              size = 3) +
    coord_equal() +
    scale_colour_manual(values = type_cols) +
    labs(
      x = expression(x[1]),
      y = expression(x[2]),
      title = "Representative measures in group C",
      colour = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  p_summary <- ggplot(summ, aes(x = group, y = estimate)) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08) +
    geom_point(aes(y = true_D_sq), shape = 4, size = 3, stroke = 1.1) +
    labs(
      x = "Group",
      y = "Estimate",
      title = expression(paste("Within-group heterogeneity for ", psi(t) == t^2))
    ) +
    theme_bw()
  
  p_eccentricity <- ggplot(ecc, aes(x = obs_index, y = ghat, colour = type)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2) +
    scale_colour_manual(values = comp_cols) +
    labs(
      x = "Observation index within group C",
      y = expression(hat(g)[i]),
      title = "Empirical eccentricities in group C",
      colour = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  p_transform <- ggplot(trans, aes(x = group, y = estimate, fill = transform)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_fill_manual(values = trans_fills) +
    labs(
      x = NULL,
      y = "Estimate",
      title = expression(paste("Effect of the transform ", psi)),
      fill = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  p_plugin <- ggplot(plugin, aes(x = inner_n, y = value, colour = series, shape = series)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = plugin_cols) +
    scale_shape_manual(values = plugin_shapes) +
    labs(
      x = "Inner sample size per measure",
      y = "Error",
      title = expression(paste("Plug-in stability for ", psi(t) == t)),
      colour = NULL,
      shape = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  list(
    locations = p_locations,
    representative = p_representative,
    summary = p_summary,
    eccentricity = p_eccentricity,
    transform = p_transform,
    plugin = p_plugin
  )
}

# -------------------------------------------------------------------------
# shared plot theme helpers for the Quarto document
# -------------------------------------------------------------------------
base_text_size   <- 16
title_text_size  <- 18
legend_text_size <- 14
axis_title_size  <- 16
axis_text_size   <- 14
strip_text_size  <- 15
text_scale <- 1.0

base_text_size   <- base_text_size   * text_scale
title_text_size  <- title_text_size  * text_scale
legend_text_size <- legend_text_size * text_scale
axis_title_size  <- axis_title_size  * text_scale
axis_text_size   <- axis_text_size   * text_scale
strip_text_size  <- strip_text_size  * text_scale

text_theme <- theme(
  text = element_text(size = base_text_size),
  axis.title = element_text(size = axis_title_size),
  axis.text = element_text(size = axis_text_size),
  legend.text = element_text(size = legend_text_size),
  legend.title = element_text(size = legend_text_size),
  strip.text = element_text(size = strip_text_size),
  plot.title = element_text(
    face = "bold",
    size = title_text_size,
    hjust = 0
  )
)

common_theme <- theme(
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.box = "horizontal",
  legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
  legend.box.margin = margin(0, 0, 0, 0),
  plot.margin = margin(3, 3, 3, 3),
  aspect.ratio = 1
)

label_panel <- function(p, label) {
  p +
    labs(title = label) +
    text_theme +
    common_theme
}

# -------------------------------------------------------------------------
# MNIST digit heterogeneity analysis
# Add this block to the end of src.R
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------------

symmetrise_distance_matrix <- function(D, tol = 1e-8) {
  stopifnot(is.matrix(D), is.numeric(D))
  if (nrow(D) != ncol(D)) stop("Distance matrix must be square.")
  if (any(!is.finite(D))) stop("Distance matrix contains non-finite values.")
  
  D <- (D + t(D)) / 2
  diag(D) <- 0
  
  if (min(D) < -tol) stop("Distance matrix has substantially negative entries.")
  D[D < 0] <- 0
  
  if (max(abs(D - t(D))) > 100 * tol) {
    stop("Distance matrix is not symmetric within tolerance.")
  }
  D
}

get_offdiag <- function(M) {
  M[upper.tri(M)]
}

load_pdmat_file <- function(file) {
  env <- new.env(parent = emptyenv())
  load(file = file, envir = env)
  if (!exists("pdmat", envir = env, inherits = FALSE)) {
    stop(sprintf("File %s does not contain an object named 'pdmat'.", file))
  }
  D <- get("pdmat", envir = env, inherits = FALSE)
  symmetrise_distance_matrix(D)
}

# -------------------------------------------------------------------------
# Kernel construction from a distance matrix
# -------------------------------------------------------------------------

make_kernel_matrix <- function(D,
                               transform = c("squared", "bounded", "identity"),
                               distance_is_squared = FALSE,
                               c0 = NULL) {
  transform <- match.arg(transform)
  D <- symmetrise_distance_matrix(D)
  
  dist_mat <- if (distance_is_squared) sqrt(pmax(D, 0)) else D
  sqdist_mat <- if (distance_is_squared) D else D^2
  
  H <- switch(
    transform,
    squared = sqdist_mat,
    identity = dist_mat,
    bounded = {
      if (is.null(c0)) {
        c0 <- stats::median(get_offdiag(dist_mat), na.rm = TRUE)
      }
      sqdist_mat / (sqdist_mat + c0^2)
    }
  )
  
  diag(H) <- 0
  attr(H, "transform") <- transform
  attr(H, "c0") <- c0
  H
}

# -------------------------------------------------------------------------
# One-sample quantities from a kernel matrix
# -------------------------------------------------------------------------

u_stat_from_kernel <- function(H, alpha = 0.05, D0 = 0) {
  stopifnot(is.matrix(H), is.numeric(H))
  n <- nrow(H)
  if (n != ncol(H)) stop("Kernel matrix must be square.")
  if (n < 3) stop("Need at least 3 observations.")
  
  U_n <- mean(H[upper.tri(H)])
  row_mean_excl <- rowSums(H) / (n - 1)
  ghat <- row_mean_excl - U_n
  sigma2_hat <- 4 * sum(ghat^2) / (n - 1)
  se <- sqrt(sigma2_hat / n)
  
  z_value <- if (se > 0) (U_n - D0) / se else NA_real_
  p_value <- if (is.na(z_value)) NA_real_ else 2 * stats::pnorm(-abs(z_value))
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  summary <- data.frame(
    n = n,
    estimate = U_n,
    sigma2_hat = sigma2_hat,
    se = se,
    ci_lower = U_n - z_crit * se,
    ci_upper = U_n + z_crit * se,
    null_value = D0,
    z_value = z_value,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
  
  ecc <- data.frame(
    obs_index = seq_len(n),
    ghat = ghat,
    rank_desc = rank(-ghat, ties.method = "first"),
    rank_asc = rank(ghat, ties.method = "first"),
    stringsAsFactors = FALSE
  )
  
  list(summary = summary, ghat = ghat, eccentricity_table = ecc)
}

# -------------------------------------------------------------------------
# Digit-level analysis
# -------------------------------------------------------------------------

analyse_one_digit <- function(D,
                              digit,
                              transform = c("squared", "bounded", "identity"),
                              distance_is_squared = FALSE,
                              alpha = 0.05,
                              D0 = 0,
                              c0 = NULL,
                              top_k = 12) {
  transform <- match.arg(transform)
  H <- make_kernel_matrix(
    D,
    transform = transform,
    distance_is_squared = distance_is_squared,
    c0 = c0
  )
  fit <- u_stat_from_kernel(H, alpha = alpha, D0 = D0)
  
  fit$summary$digit <- as.character(digit)
  fit$summary$transform <- transform
  fit$summary$c0 <- if (is.null(attr(H, "c0"))) NA_real_ else attr(H, "c0")
  
  ecc <- fit$eccentricity_table
  ecc$digit <- as.character(digit)
  ecc$transform <- transform
  
  top_ecc <- ecc[order(-ecc$ghat, ecc$obs_index), , drop = FALSE]
  top_ecc <- head(top_ecc, min(top_k, nrow(top_ecc)))
  top_ecc$which <- "top_eccentric"
  
  central <- ecc[order(ecc$ghat, ecc$obs_index), , drop = FALSE]
  central <- head(central, min(top_k, nrow(central)))
  central$which <- "most_central"
  
  list(
    summary = fit$summary,
    eccentricity_table = ecc,
    top_eccentric = top_ecc,
    most_central = central,
    kernel_matrix = H
  )
}

compare_two_fits <- function(fitA, fitB, digitA, digitB, alpha = 0.05, diff0 = 0) {
  sA <- fitA$summary
  sB <- fitB$summary
  
  diff_est <- sA$estimate - sB$estimate
  se_diff <- sqrt(sA$sigma2_hat / sA$n + sB$sigma2_hat / sB$n)
  z_value <- if (se_diff > 0) (diff_est - diff0) / se_diff else NA_real_
  p_value <- if (is.na(z_value)) NA_real_ else 2 * stats::pnorm(-abs(z_value))
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  data.frame(
    digit_A = as.character(digitA),
    digit_B = as.character(digitB),
    estimate_A = sA$estimate,
    estimate_B = sB$estimate,
    estimate_diff = diff_est,
    se_diff = se_diff,
    ci_lower = diff_est - z_crit * se_diff,
    ci_upper = diff_est + z_crit * se_diff,
    null_value = diff0,
    z_value = z_value,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

analyse_mnist_digits <- function(data_dir = ".",
                                 digits = 0:9,
                                 distance_is_squared = FALSE,
                                 main_transform = c("squared", "bounded", "identity"),
                                 do_bounded_sensitivity = TRUE,
                                 bounded_c0 = NULL,
                                 bounded_c0_rule = c("global_median", "per_digit_median"),
                                 alpha = 0.05,
                                 D0 = 0,
                                 top_k = 12) {
  main_transform <- match.arg(main_transform)
  bounded_c0_rule <- match.arg(bounded_c0_rule)
  
  D_list <- list()
  for (d in digits) {
    file <- file.path(data_dir, sprintf("digit_%s.RData", d))
    if (!file.exists(file)) stop(sprintf("Missing file: %s", file))
    D_list[[as.character(d)]] <- load_pdmat_file(file)
  }
  
  if (do_bounded_sensitivity && is.null(bounded_c0) && bounded_c0_rule == "global_median") {
    pooled_distances <- unlist(lapply(D_list, function(D) {
      dist_mat <- if (distance_is_squared) sqrt(pmax(D, 0)) else D
      get_offdiag(dist_mat)
    }), use.names = FALSE)
    bounded_c0 <- stats::median(pooled_distances, na.rm = TRUE)
  }
  
  main_fits <- list()
  summary_main <- list()
  ecc_main <- list()
  top_main <- list()
  central_main <- list()
  
  for (nm in names(D_list)) {
    c0_use <- NULL
    if (main_transform == "bounded") {
      c0_use <- if (!is.null(bounded_c0)) bounded_c0 else NULL
      if (is.null(c0_use) && bounded_c0_rule == "per_digit_median") {
        dist_mat <- if (distance_is_squared) sqrt(pmax(D_list[[nm]], 0)) else D_list[[nm]]
        c0_use <- stats::median(get_offdiag(dist_mat), na.rm = TRUE)
      }
    }
    
    fit <- analyse_one_digit(
      D = D_list[[nm]],
      digit = nm,
      transform = main_transform,
      distance_is_squared = distance_is_squared,
      alpha = alpha,
      D0 = D0,
      c0 = c0_use,
      top_k = top_k
    )
    
    main_fits[[nm]] <- fit
    summary_main[[nm]] <- fit$summary
    ecc_main[[nm]] <- fit$eccentricity_table
    top_main[[nm]] <- fit$top_eccentric
    central_main[[nm]] <- fit$most_central
  }
  
  summary_main_df <- do.call(rbind, summary_main)
  rownames(summary_main_df) <- NULL
  summary_main_df$digit <- factor(summary_main_df$digit, levels = as.character(digits))
  summary_main_df$rank_desc <- rank(-summary_main_df$estimate, ties.method = "first")
  summary_main_df$rank_asc <- rank(summary_main_df$estimate, ties.method = "first")
  
  ecc_main_df <- do.call(rbind, ecc_main)
  rownames(ecc_main_df) <- NULL
  ecc_main_df$digit <- factor(ecc_main_df$digit, levels = as.character(digits))
  
  top_main_df <- do.call(rbind, top_main)
  rownames(top_main_df) <- NULL
  top_main_df$digit <- factor(top_main_df$digit, levels = as.character(digits))
  
  central_main_df <- do.call(rbind, central_main)
  rownames(central_main_df) <- NULL
  central_main_df$digit <- factor(central_main_df$digit, levels = as.character(digits))
  
  pairwise_list <- list()
  idx <- 1L
  for (i in seq_along(digits)) {
    for (j in seq_along(digits)) {
      if (i < j) {
        di <- as.character(digits[i])
        dj <- as.character(digits[j])
        pairwise_list[[idx]] <- compare_two_fits(
          fitA = main_fits[[di]],
          fitB = main_fits[[dj]],
          digitA = di,
          digitB = dj,
          alpha = alpha,
          diff0 = 0
        )
        idx <- idx + 1L
      }
    }
  }
  pairwise_df <- do.call(rbind, pairwise_list)
  rownames(pairwise_df) <- NULL
  pairwise_df$p_value_bh <- stats::p.adjust(pairwise_df$p_value, method = "BH")
  
  bounded_summary_df <- NULL
  transform_compare_df <- NULL
  bounded_fits <- NULL
  
  if (do_bounded_sensitivity && main_transform != "bounded") {
    bounded_fits <- list()
    bounded_summary <- list()
    
    for (nm in names(D_list)) {
      c0_use <- bounded_c0
      if (is.null(c0_use) && bounded_c0_rule == "per_digit_median") {
        dist_mat <- if (distance_is_squared) sqrt(pmax(D_list[[nm]], 0)) else D_list[[nm]]
        c0_use <- stats::median(get_offdiag(dist_mat), na.rm = TRUE)
      }
      
      fit_b <- analyse_one_digit(
        D = D_list[[nm]],
        digit = nm,
        transform = "bounded",
        distance_is_squared = distance_is_squared,
        alpha = alpha,
        D0 = D0,
        c0 = c0_use,
        top_k = top_k
      )
      bounded_fits[[nm]] <- fit_b
      bounded_summary[[nm]] <- fit_b$summary
    }
    
    bounded_summary_df <- do.call(rbind, bounded_summary)
    rownames(bounded_summary_df) <- NULL
    bounded_summary_df$digit <- factor(bounded_summary_df$digit, levels = as.character(digits))
    
    transform_compare_df <- merge(
      summary_main_df[, c("digit", "estimate")],
      bounded_summary_df[, c("digit", "estimate")],
      by = "digit",
      suffixes = c("_main", "_bounded")
    )
    transform_compare_df$ratio_bounded_to_main <- with(
      transform_compare_df,
      estimate_bounded / estimate_main
    )
  }
  
  list(
    summary = summary_main_df,
    pairwise = pairwise_df,
    eccentricity = ecc_main_df,
    top_eccentric = top_main_df,
    most_central = central_main_df,
    bounded_summary = bounded_summary_df,
    transform_compare = transform_compare_df,
    fits = main_fits,
    bounded_fits = bounded_fits,
    matrices = D_list,
    settings = list(
      distance_is_squared = distance_is_squared,
      main_transform = main_transform,
      alpha = alpha,
      D0 = D0,
      bounded_c0 = bounded_c0,
      bounded_c0_rule = bounded_c0_rule,
      top_k = top_k
    )
  )
}

# -------------------------------------------------------------------------
# Optional output helpers
# -------------------------------------------------------------------------

plot_digit_summary <- function(summary_df,
                               title = expression(paste("Estimated digit-wise heterogeneity for ", psi(t) == t^2))) {
  summary_df$digit <- factor(summary_df$digit, levels = as.character(0:9))
  
  ggplot(summary_df, aes(x = digit, y = estimate)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
    geom_point(size = 2.2) +
    labs(x = "Digit", y = "Estimate", title = title) +
    theme_bw()
}

plot_transform_comparison <- function(transform_compare_df,
                                      title = "Squared versus bounded transform") {
  long_df <- rbind(
    data.frame(
      digit = transform_compare_df$digit,
      transform = "squared",
      estimate = transform_compare_df$estimate_main
    ),
    data.frame(
      digit = transform_compare_df$digit,
      transform = "bounded",
      estimate = transform_compare_df$estimate_bounded
    )
  )
  long_df$digit <- factor(long_df$digit, levels = as.character(0:9))
  
  ggplot(long_df, aes(x = digit, y = estimate, fill = transform)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    labs(x = "Digit", y = "Estimate", title = title, fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

plot_top_eccentric_values <- function(top_eccentric_df,
                                      top_n_per_digit = 8,
                                      title = "Largest empirical eccentricities by digit") {
  keep <- do.call(rbind, lapply(split(top_eccentric_df, top_eccentric_df$digit), function(df) {
    df <- df[order(-df$ghat), , drop = FALSE]
    head(df, min(top_n_per_digit, nrow(df)))
  }))
  rownames(keep) <- NULL
  keep$digit <- factor(keep$digit, levels = as.character(0:9))
  keep$rank_within_digit <- ave(-keep$ghat, keep$digit, FUN = seq_along)
  
  ggplot(keep, aes(x = rank_within_digit, y = ghat, group = digit)) +
    geom_line() +
    geom_point(size = 1.7) +
    facet_wrap(~ digit, scales = "free_y") +
    labs(x = "Rank within digit", y = expression(hat(g)[i]), title = title) +
    theme_bw()
}

extract_mnist_images_by_index <- function(images_one_digit, indices) {
  if (length(dim(images_one_digit)) == 3) {
    return(lapply(indices, function(i) images_one_digit[, , i]))
  }
  if (is.matrix(images_one_digit) && ncol(images_one_digit) == 28 * 28) {
    return(lapply(indices, function(i) matrix(images_one_digit[i, ], nrow = 28, byrow = TRUE)))
  }
  stop("images_one_digit must be either a 28 x 28 x n array or an n x 784 matrix.")
}