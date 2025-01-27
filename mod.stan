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
    matrix [N_lin, 2] site_gyrA_presence_absence;
    matrix [N_lin, 4] site_penA_presence_absence;
    matrix [N_lin, 1] site_pbp_presence_absence;
    matrix [N_lin, 5] site_mtr_presence_absence;
    matrix [N_lin, 1] site_rRNA23S_presence_absence;

    array[N_lin] int<lower=0> site_gyrA_re_idx;
    matrix[N_lin, 2] site_gyrA_Z;

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
    vector[2] site_gyrA_a_ncp;
    vector[4] site_penA_a_ncp;
    vector[1] site_pbp_a_ncp;
    vector[5] site_mtr_a_ncp;
    vector[1] site_rRNA23S_a_ncp;

    //SITE EFFECT PARAMETERS
    matrix [2, 1] site_gyrA_b_ncp;
    matrix [4, 1] site_penA_b_ncp;
    matrix [1, 1] site_pbp_b_ncp;
    matrix [1, 1] site_rRNA23S_b_ncp;
    matrix [5, 3] site_mtr_b_ncp;
    vector [5] site_mtr_b_eff;
    real <lower=0, upper=1> site_mtr_pool;

    array[2, 3, 1] real site_gyrA_re_b_raw;
    matrix[2, 3] site_gyrA_re_a_raw;
    real <lower=0> site_gyrA_re_sd_b;
    real <lower=0> site_gyrA_re_sd_a;

}

transformed parameters {

    real <lower=0> offset_intercept_sd = sc_factor * offset_intercept_sd_sc;
    real <lower=0> tau_trend = sc_factor * tau_trend_sc;
    real <lower=0> group_sd = sc_factor * group_sd_sc;
    real <lower=0> residual_sd = sc_factor * residual_sd_sc / sqrt(grid_width);

    //SITE BS
    matrix [5, 3] site_mtr_b = site_mtr_pool * rep_matrix(site_mtr_b_eff * 0.25, 3) + (1.0 - site_mtr_pool) * site_mtr_b_ncp * 0.25;
    matrix [2, 1] site_gyrA_b = 0.25 * site_gyrA_b_ncp;
    matrix [4, 1] site_penA_b = 0.25 * site_penA_b_ncp;
    matrix [1, 1] site_pbp_b = 0.25 * site_pbp_b_ncp;
    matrix [1, 1] site_rRNA23S_b = 0.25 * site_rRNA23S_b_ncp;
    
    //SITE INTERCEPTS
    vector [2] site_gyrA_a = 0.25 * site_gyrA_a_ncp;
    vector [4] site_penA_a = 0.25 * site_penA_a_ncp;
    vector [1] site_pbp_a = 0.25 * site_pbp_a_ncp;
    vector [5] site_mtr_a = 0.25 * site_mtr_a_ncp;
    vector [1] site_rRNA23S_a = 0.25 * site_rRNA23S_a_ncp;

    array[2] matrix[3, 1] site_gyrA_re_b = scale_re_b(site_gyrA_re_b_raw, site_gyrA_re_sd_b * 0.25);
    matrix[2, 3] site_gyrA_re_a = site_gyrA_re_sd_a * 0.25 * site_gyrA_re_a_raw;

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
    vector[N_lin] site_gyrA_lineage_a = site_gyrA_presence_absence * site_gyrA_a;
    matrix[N_lin, 1] site_gyrA_lineage_b = site_gyrA_presence_absence * site_gyrA_b;
    vector[N_lin] site_penA_lineage_a = site_penA_presence_absence * site_penA_a;
    matrix[N_lin, 1] site_penA_lineage_b = site_penA_presence_absence * site_penA_b;
    vector[N_lin] site_pbp_lineage_a = site_pbp_presence_absence * site_pbp_a;
    matrix[N_lin, 1] site_pbp_lineage_b = site_pbp_presence_absence * site_pbp_b;
    vector[N_lin] site_mtr_lineage_a = site_mtr_presence_absence * site_mtr_a;
    matrix[N_lin, 3] site_mtr_lineage_b = site_mtr_presence_absence * site_mtr_b;
    vector[N_lin] site_rRNA23S_lineage_a = site_rRNA23S_presence_absence * site_rRNA23S_a;
    matrix[N_lin, 1] site_rRNA23S_lineage_b = site_rRNA23S_presence_absence * site_rRNA23S_b;

    matrix [N_lin, 1] re_gyrA_lin_b = lin_re_b(site_gyrA_re_b, site_gyrA_Z, site_gyrA_re_idx);
    vector [N_lin] re_gyrA_lin_a = lin_re_a(site_gyrA_re_a, site_gyrA_Z, site_gyrA_re_idx);

    //ALLELE RT EFFECTS
    lin_trends += rep_matrix(site_gyrA_lineage_a, N_grid);
    lin_trends += rep_matrix(site_penA_lineage_a, N_grid);
    lin_trends += rep_matrix(site_pbp_lineage_a, N_grid);
    lin_trends += rep_matrix(site_mtr_lineage_a, N_grid);
    lin_trends += rep_matrix(site_rRNA23S_lineage_a, N_grid);
    lineage_coeffs[:, 1] += site_gyrA_lineage_b[:, 1];
    lineage_coeffs[:, 2] += site_penA_lineage_b[:, 1];
    lineage_coeffs[:, 2] += site_pbp_lineage_b[:, 1];
    lineage_coeffs[:, 1] += site_mtr_lineage_b[:, 1];
    lineage_coeffs[:, 2] += site_mtr_lineage_b[:, 2];
    lineage_coeffs[:, 3] += site_mtr_lineage_b[:, 3];
    lineage_coeffs[:, 3] += site_rRNA23S_lineage_b[:, 1];

    //RANDOM RT EFFECTS
    lin_trends += rep_matrix(re_gyrA_lin_a, N_grid);
    lineage_coeffs[:, 1] += re_gyrA_lin_b[:, 1];

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
    site_gyrA_a_ncp ~ normal(0, 4);
    site_penA_a_ncp ~ normal(0, 4);
    site_pbp_a_ncp ~ normal(0, 4);
    site_mtr_a_ncp ~ normal(0, 4);
    site_rRNA23S_a_ncp ~ normal(0, 4);
    to_array_1d(site_gyrA_b_ncp) ~ normal(0, 4);
    to_array_1d(site_penA_b_ncp) ~ normal(0, 4);
    to_array_1d(site_pbp_b_ncp) ~ normal(0, 4);
    to_array_1d(site_rRNA23S_b_ncp) ~ normal(0, 4);
    site_mtr_b_eff ~ normal(0, 4);
    to_array_1d(site_mtr_b_ncp) ~ normal(0, 4);

    //RANDOM EFFECT PRIORS
    to_array_1d(site_gyrA_re_a_raw) ~ normal(0, 1);
    to_array_1d(site_gyrA_re_b_raw) ~ normal(0, 1);
    site_gyrA_re_sd_a ~ normal(0, 4);
    site_gyrA_re_sd_b ~ normal(0, 4);

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
    lineage_intercepts += site_gyrA_lineage_a;
    lineage_intercepts += site_penA_lineage_a;
    lineage_intercepts += site_pbp_lineage_a;
    lineage_intercepts += site_mtr_lineage_a;
    lineage_intercepts += site_rRNA23S_lineage_a;
}