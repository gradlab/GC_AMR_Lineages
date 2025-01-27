functions {
    real coal_lik(vector f, vector aik, vector laik, vector wt, array[] int ce, array[] int g_idx) {
        real out = 0.0;
        out += -sum(aik .* wt ./ f[g_idx]);
        out += -sum(log(f[g_idx[ce]]));
        return out;
    }

    array[] matrix scale_re_b(array[,,] real re_raw, real sd_re) {
        int nrow = dims(re_raw)[2];
        int ncol = dims(re_raw)[3];
        int N = dims(re_raw)[1];

        array[N] matrix[nrow, ncol] out;

        for (i in 1:N) {
            out[i] = sd_re * to_matrix(re_raw[i,:,:]);
        }

        return out;
    }

    matrix lin_re_b(array[] matrix re_b, matrix re_Z, array[] int grp_idx) {
        int N_motif = size(re_b);
        int N_lin = size(grp_idx);
        int N_covs = cols(re_b[1]);

        matrix[N_lin, N_covs] out = rep_matrix(0.0, N_lin, N_covs);

        for (i in 1:N_lin) {
            for (j in 1:N_motif) {
                int re_i = grp_idx[i];
                if (re_i != 0) {
                    out[i, :] += re_Z[i, j] * re_b[j][re_i, :];
                }
            }
        }

        return out;
    }

    vector lin_re_a(matrix re_a, matrix re_Z, array[] int grp_idx) {
        int N_lvls = cols(re_a);
        int N_motif = rows(re_a);
        int N_lin = size(grp_idx);
        vector[N_lin] out = rep_vector(0.0, N_lin);
        
        for (i in 1:N_lin) {
            for (j in 1:N_motif) {
                int re_i = grp_idx[i];
                if (re_i != 0) {
                    out[i] += re_Z[i, j] * re_a[j, re_i];
                }
            }
        }
        return out;
    }
    
    vector rw_ncp(vector ft, real tau, real grid_width) {
        int N = size(ft);
        vector[N] temp;
        temp[1] = 0.0;
        temp[2:N] = ft[2:N] * tau * sqrt(grid_width);
        return cumulative_sum(temp);
    }

    row_vector rw_ncp(row_vector ft, real tau, real grid_width) {
        int N = size(ft);
        row_vector[N] temp;
        temp[1] = 0.0;
        temp[2:N] = ft[2:N] * tau * sqrt(grid_width);
        return cumulative_sum(temp);
    }

    vector mean_traj(vector rt, real N0, real width) {
        return(N0 * exp(cumulative_sum(-width * rt)));
    }

    matrix mean_traj(matrix rts, vector N0s, real width)
    {
        int nrow = rows(rts);
        int ncol = cols(rts);
        matrix[nrow, ncol] out;

        for (i in 1:nrow)
        {
            out[i, :] = N0s[i] * exp(cumulative_sum(-width * rts[i, :]));
        }

        return out;
    }

    matrix make_rt_intercept_vcv_ld(int N_lin, real group_sd, real mu_sd) {
        matrix[N_lin+1, N_lin+1] rt_intercept_vcv = rep_matrix(mu_sd * mu_sd, N_lin+1, N_lin+1);
        vector[N_lin+1] diag = rep_vector(1e-8, N_lin+1);
        diag[1:N_lin] += group_sd * group_sd;
        rt_intercept_vcv = add_diag(rt_intercept_vcv, diag);
        
        return cholesky_decompose(rt_intercept_vcv);
    }

    matrix vec_arr_to_matrix(array[] vector vec_arr, int n, int m){
        matrix[n, m] out;
        for (i in 1:n)
        {
            out[i, :] = to_row_vector(vec_arr[i]);
        }
        return out;
    }
}



data {
    int<lower=0> N_tot_obs;
    int<lower=0> N_tot_coal;
    int<lower=0> N_grid;
    int<lower=0> N_lin; 
    array[N_lin] int<lower=0> lin_sizes;
    array[N_lin] int<lower=0> lin_coal_counts;

    vector[N_tot_obs] wt;
    vector[N_tot_obs] aik;
    array[N_tot_coal] int<lower=0> ce;
    array[N_tot_obs] int<lower=0> grid_idx;
    real<lower=0> grid_width;
    
    int<lower=0> K_covar;
    matrix[K_covar, N_grid] covars;

    //ALLELE-COVARIATE INTERACTION MATRICES
###```DATA_BLOCK```

###```RE_GRP_BLOCK```

}

transformed data {
    vector[N_tot_obs] laik = log(aik);  
    real sc_factor = 1.0 / 8.0;
}

parameters {    
    vector[N_lin+1] rt_intercepts_ncp;
    vector[N_grid] global_trend_ncp;
    
    vector[K_covar] beta_global;
    
    vector[N_lin] lin_offsets_ncp;
    real offset_mu;
    
    real <lower=0> offset_intercept_sd_sc;
    real <lower=0> tau_trend_sc;
    real <lower=0> group_sd_sc;
    
    real <lower=0> residual_sd_sc;
    matrix[N_lin, N_grid] f_residual;

    //SITE INTERCEPT PARAMETERS
###```SITE_INTERCEPT_BLOCK```

    //SITE EFFECT PARAMETERS
###```B_RAW_BLOCK```

###```RE_BLOCK```

}

transformed parameters {

    real <lower=0> offset_intercept_sd = sc_factor * offset_intercept_sd_sc;
    real <lower=0> tau_trend = sc_factor * tau_trend_sc;
    real <lower=0> group_sd = sc_factor * group_sd_sc;
    real <lower=0> residual_sd = sc_factor * residual_sd_sc / sqrt(grid_width);

    //SITE BS
###```B_TRANSF_BLOCK```
    
    //SITE INTERCEPTS
###```A_TRANSF_BLOCK```

###```RE_TRANSF_BLOCK```

    matrix[N_lin, N_grid] lin_trends = rep_matrix(0.0, N_lin, N_grid);

    vector[N_grid] global_trend;
    vector[N_grid] global_effect;
    vector[N_lin] lin_rt_intercepts;
    real rt_mu;
    
    {
        global_effect = to_vector(to_row_vector(beta_global * 0.5) * covars);

        matrix[N_lin + 1, N_lin + 1] rt_intercept_vcv_ld = make_rt_intercept_vcv_ld(N_lin, group_sd, 1.0);
        vector[N_lin + 1] rt_intercepts = rt_intercept_vcv_ld * rt_intercepts_ncp;
        
        rt_mu = rt_intercepts[N_lin + 1];
        lin_rt_intercepts = rt_intercepts[1:N_lin] - rt_mu;
        global_trend = rw_ncp(global_trend_ncp, tau_trend, grid_width) + rt_mu;

        for (i in 1:N_lin)
        {
            lin_trends[i,:] += lin_rt_intercepts[i];
        }
    }
 
    vector[N_lin] lin_offsets = exp(lin_offsets_ncp * offset_intercept_sd + offset_mu + 4.0);
    matrix[N_lin, N_grid] lin_traj = rep_matrix(0.0, N_lin, N_grid);

    matrix[N_lin, K_covar] lineage_coeffs = rep_matrix(0.0, N_lin, K_covar);

    //ALLELE BS
###```TPARAMS_BLOCK```

###```RE_LIN_BLOCK```

    //ALLELE RT EFFECTS
###```EFFECTS_BLOCK```

    //RANDOM RT EFFECTS
###```RE_EFF_BLOCK```

    lin_trends += lineage_coeffs * covars;

    matrix [N_lin, N_grid] lin_rt_resi = residual_sd * f_residual;

    //COMPUTE TRAJECTORIES
    for (i in 1:N_lin) {
        lin_rt_resi[i, :] += lin_trends[i,:] + to_row_vector(global_trend + global_effect);
    }

    lin_traj = mean_traj(lin_rt_resi, lin_offsets, grid_width);
}

model {

    beta_global ~ normal(0, 2);
    rt_intercepts_ncp ~ normal(0, 1);
    lin_offsets_ncp ~ normal(0, 1);
    
    global_trend_ncp ~ normal(0, 1);
    tau_trend ~ normal(0, 1);
    offset_mu ~ normal(0, 4);
    offset_intercept_sd ~ normal(0, 1);

    group_sd ~ normal(0, 1);

    //RESIDUAL PRIORS
    to_array_1d(f_residual) ~ normal(0, 1);
    residual_sd ~ normal(0,1);

    //ALLELE EFFECT PRIORS
###```MODEL_BLOCK```

    //RANDOM EFFECT PRIORS
###```RE_MODEL_BLOCK```

    int sz_idx = 1;
    int ce_idx = 1;

    for (i in 1:N_lin) {
        int lsz = lin_sizes[i];
        int lcc = lin_coal_counts[i];
        vector[lsz] lin_wt = segment(wt, sz_idx, lsz);
        vector[lsz] lin_aik = segment(aik, sz_idx, lsz);
        vector[lsz] lin_laik = segment(laik, sz_idx, lsz);
        array[lsz] int lin_grid_idx = segment(grid_idx, sz_idx, lsz);
        array[lcc] int lin_ce = segment(ce, ce_idx, lcc);

        target += coal_lik(to_vector(lin_traj[i,:]), lin_aik, lin_laik, lin_wt, lin_ce, lin_grid_idx);

        sz_idx += lsz;
        ce_idx += lcc;
    }
}

generated quantities {
    matrix [N_lin, N_grid] lin_rt_noresi;
    matrix [N_lin, N_grid] lin_traj_noresi;

    vector [N_lin] lineage_intercepts = lin_rt_intercepts;

    for (i in 1:N_lin) {
        lin_rt_noresi[i, :] = lin_trends[i,:] + to_row_vector(global_trend + global_effect);
        lin_traj_noresi[i, :] = to_row_vector(mean_traj(to_vector(lin_rt_noresi[i, :]), lin_offsets[i], grid_width));
    }
        //LINEAGE EFFECTS AND INTERCEPTS
###```GENERATED_BLOCK```
}