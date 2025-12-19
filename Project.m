clc; clear;
alpha = (-5:20)';
CL_2412 = [-0.56286, -0.36946, -0.18122, -0.075955, 0.10661, 0.31453, 0.53035, 0.71075, 0.90488, 1.1854, 1.1897, 1.4191, 1.6484, 1.8747, 2.0811, 2.3373, 2.5595, 2.821, 3.1487, 3.1738, 3.5495, 3.6152, 4.1389, 4.4859, 6.0614, 4.9339]';
CL_2414 = [-0.57526, -0.38848, -0.19561, -0.082528, 0.12287, 0.30777, 0.48876, 0.67904, 0.83113, 1.0554, 1.2107, 1.4091, 1.624, 1.8471, 2.1323, 2.3274, 2.585, 2.8738, 3.1244, 3.441, 3.6676, 4.2174, 4.1574, 4.5851, 5.0823, 6.1681]';
CL_2415 = [-0.68794, -0.48322, -0.2876, -0.092863, 0.10895, 0.33634, 0.48311, 0.8445, 0.89864, 1.1284, 1.3402, 1.6187, 1.8308, 1.8623, 2.1178, 2.5307, 2.5667, 2.8194, 3.1779, 3.4854, 3.9596, 4.1032, 4.0981, 4.617, 4.4271, 4.906]';
y = CL_2412; % selected airfoil for cur
m = length(alpha);

%% Givens Rotations Analysis
% Error metrics
rmse    = @(y, yhat) sqrt(mean((y - yhat).^2));
r_coeff = @(y, yhat) 1 - sum((y-yhat).^2) / sum((y-mean(y)).^2);

% Linear fit using Givens Rottion function
X_linear = [ones(size(alpha)) alpha];
tic; % start of timer 
[R_linear, QTy_linear] = givens(X_linear, y);
linear_coeff   = R_linear \ QTy_linear;
yhat_linear    = X_linear * linear_coeff;
runtime_linear = toc; % end of timer

% Metrics
rmse_linear    = rmse(y, yhat_linear);
r_coeff_linear = r_coeff(y, yhat_linear);
cond_linear    = cond(X_linear);

% Flop counts
c = 2; % 2 columns in linear model
flops_linear = 6*m*c - 3*c*(c+1);

% Polynomial fits: degrees 2–6 using Givens
degrees_fit = 2:6;
rmse_poly    = zeros(length(degrees_fit), 1);
r2_poly      = zeros(length(degrees_fit), 1);
cond_poly    = zeros(length(degrees_fit), 1);
runtime_poly = zeros(length(degrees_fit), 1);
flops_poly   = zeros(length(degrees_fit), 1);

for i = 1:length(degrees_fit)
    n = degrees_fit(i);
    
    % Vandermonde matrix
    X_poly = zeros(length(alpha), n+1);
    for k = 0:n
        X_poly(:, k+1) = alpha.^k;
    end
    
    % QR via Givens implementation
    tic;
    [R_poly, QTy_poly] = givens(X_poly, y);
    poly_coeff       = R_poly \ QTy_poly;
    yhat_poly        = X_poly * poly_coeff;
    runtime_poly(i)  = toc;
    
    % Metrics
    rmse_poly(i) = rmse(y, yhat_poly);
    r2_poly(i)   = r_coeff(y, yhat_poly);
    cond_poly(i) = cond(X_poly);
    
    % Flop counts
    c = n + 1;
    flops_poly(i) = 6*m*c - 3*c*(c+1);
    
    % Plot data, linear, and polynomial fit
    figure(i);
    plot(alpha, y, 'o', alpha, yhat_linear, '-', alpha, yhat_poly, '--');
    grid on;
    title(sprintf('Polynomial Degree %d Comparison', n));
    legend('Data', 'Linear fit', sprintf('Poly-%d fit', n), 'Location', 'best');
    xlabel('\alpha (degrees)');
    ylabel('C_L');
end

% Summary tables
linear_table = table(rmse_linear, r_coeff_linear, cond_linear, runtime_linear, flops_linear);
poly_table   = table(degrees_fit', rmse_poly, r2_poly, cond_poly, runtime_poly, flops_poly);
disp(linear_table)
disp(poly_table)

%% Leave one out cross validation
% set CV storages
degrees_CV = 1:24;
ndeg = length(degrees_CV);
yhat_CV  = zeros(m, ndeg);
resid_CV = zeros(m, ndeg);
rmse_CV  = zeros(ndeg, 1);

% Coefficients
maxc = max(degrees_CV) + 1;
coeff_CV = NaN(maxc, m, ndeg);

for d = 1:length(degrees_CV)
    n = degrees_CV(d);
    c = n + 1;

    for j = 1:m

        % Build new data sets leaving out one point
        alpha_new = alpha;
        alpha_new(j) = [];
        y_new = y;
        y_new(j) = [];

        % Vandermonde matrix on new datat set
        X_new = zeros(m-1, c);
        for k = 0:n
            X_new(:, k+1) = alpha_new.^k;
        end

        % Fit with givens function
        [R_new, QTy_new] = givens(X_new, y_new);
        coeff_new = R_new \ QTy_new;
        coeff_CV(1:c, j, d) = coeff_new;

        % Prediction of left out point
        x_pre = zeros(1, c);
        for k = 0:n
            x_pre(k+1) = alpha(j)^k;
        end

        % Find RMSE for this degree
        yhat_CV(j, d) = x_pre * coeff_new;
        resid_CV(j, d) = y(j) - yhat_CV(j, d);
    end
    rmse_CV(d) = sqrt(mean(resid_CV(:, d).^2));
end

% Find coefficient variance across fits
coeff_var = NaN(maxc, ndeg);
coeff_mean = NaN(maxc, ndeg);
for d = 1:ndeg
    n = degrees_CV(d);
    c = n + 1;
    th_coeffs = squeeze(coeff_CV(1:c, :, d))';
    coeff_var(1:c, d) = var(th_coeffs, 0, 1)';
    % coeff_mean(1:c, d) = mean(th_coeffs, 1)';
end

% Results tables
CV_table = table(degrees_CV', rmse_CV);
disp(CV_table)
var_table = table(coeff_var);
disp(var_table)

%% Stability Comparisons: Householder vs Cholesky vs Givens
degrees_comp = (2:25); % change range to see exactly where each method breaks down
num_deg = length(degrees_comp);

% set RMSE for each method (matlab qr is for reference)
rmse_QR = zeros(num_deg, 1);
rmse_HH = zeros(num_deg, 1);
rmse_CH = zeros(num_deg, 1);
rmse_GR = zeros(num_deg, 1);
% set relative error for the three compared methods
err_HH = zeros(num_deg, 1);
err_CH = zeros(num_deg, 1);
err_GR = zeros(num_deg, 1);
% cond(X) for each degree
condX_comp = zeros(num_deg, 1);

for idx = 1:num_deg
    n = degrees_comp(idx);

    % Vandermode matrix
    X_comp = zeros(m, n+1);
    for k = 0:n
        X_comp(:, k+1) = alpha.^k;
    end
    condX_comp(idx) = cond(X_comp);

    % Matlab QR reference
    [Q_QR, R_QR] = qr(X_comp, 0);
    beta_QR      = R_QR \ (Q_QR' * y);
    yhat_QR      = X_comp * beta_QR;
    rmse_QR(idx) = rmse(y, yhat_QR);

    % Householder Reflections via qrfact from textbook
    [Q_HH, R_HH] = qrfact(X_comp);
    beta_HH      = R_HH \ (Q_HH' * y);
    yhat_HH      = X_comp * beta_HH;
    rmse_HH(idx) = rmse(y, yhat_HH);
    err_HH(idx)  = norm(beta_HH - beta_QR) / norm(beta_QR);

    % Cholesky (normal equations
    ATA = X_comp' * X_comp;
    ATy = X_comp' * y;
    try
        R_CH = chol(ATA);
        beta_CH      = R_CH \ (R_CH' \ ATy);
        yhat_CH      = X_comp * beta_CH;
        rmse_CH(idx) = rmse(y, yhat_CH);
        err_CH(idx)  = norm(beta_CH - beta_QR) / norm(beta_QR);
    catch
        rmse_CH(idx) = NaN;
        err_CH(idx) = NaN;
    end

    % Givens rotation QR
    [R_GR, QTy_GR] = givens(X_comp, y);
    beta_GR        = R_GR \ QTy_GR;
    yhat_GR        = X_comp * beta_GR;
    rmse_GR(idx)   = rmse(y, yhat_GR);
    err_GR(idx)    = norm(beta_GR - beta_QR) / norm(beta_QR);
end

% Table summarizing where each method derives from Matlab QR
method_table = table(degrees_comp', condX_comp, rmse_QR, rmse_HH, rmse_CH, rmse_GR, err_HH, err_CH, err_GR);
disp(method_table)

% Plots for report
figure;
semilogy(degrees_comp, condX_comp, '-o');
grid on;
xlabel('Polynomial degree');
ylabel('cond(X)');
title('Condition number vs degree');

figure;
semilogy(degrees_comp, err_HH, '-o', degrees_comp, err_CH, '-s', degrees_comp, err_GR, '-^');
grid on;
xlabel('Polynomial degree');
ylabel('Relative coefficient error vs MATLAB qr');
legend('Householder','Cholesky','Givens','Location','best');
title('Breakdown of methods relative to MATLAB qr');