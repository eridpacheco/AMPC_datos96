function [S, stats] = run_adaptivity_figures_datos96_unified(matFile, varargin)
%RUN_ADAPTIVITY_FIGURES_DATOS96_UNIFIED Generate adaptivity figures for datos_96.
%
% Mathematical definitions (for Methods section):
% Given the identified parameters at each sample k, define the parameter
% vector as:
%   theta(k) = [A11(k), A21(k), A12(k), A22(k), B1(k), B2(k)].
% The parameter increment is:
%   Delta theta(k) = theta(k) - theta(k-1),  k >= 2,
% and Delta theta(1) = NaN to preserve alignment.
% The magnitude of adaptation is the Euclidean norm:
%   ||Delta theta(k)||_2 = sqrt( sum_i (Delta theta_i(k))^2 ).
% Daily aggregation uses calendar-day bins defined by
%   day(k) = dateshift(t(k), 'start', 'day'),
% and the daily statistic is the median (or mean if configured).
%
% Usage:
%   [S, stats] = run_adaptivity_figures_datos96('datos_96_act.mat', ...
%       'SaveDir', 'figs', 'DailyAgg', 'median', ...
%       'StartDate', '2024-01-12 00:00:00', ...
%       'EndDate', '2024-02-05 23:59:59');

p = inputParser;
p.addRequired('matFile', @(x) ischar(x) || isstring(x));
p.addParameter('SaveDir', 'figs', @(x) ischar(x) || isstring(x));
p.addParameter('DailyAgg', 'median', @(x) ischar(x) || isstring(x));
p.addParameter('ShadeIrrigation', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('UseSemilogDtheta', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('StartDate', [], @(x) isempty(x) || ischar(x) || isstring(x) || isnumeric(x) || isdatetime(x));
p.addParameter('EndDate', [], @(x) isempty(x) || ischar(x) || isstring(x) || isnumeric(x) || isdatetime(x));
p.parse(matFile, varargin{:});

saveDir = char(p.Results.SaveDir);
dailyAgg = lower(char(p.Results.DailyAgg));
shadeIrrigation = logical(p.Results.ShadeIrrigation);
useSemilogDtheta = logical(p.Results.UseSemilogDtheta);
startDateInput = p.Results.StartDate;
endDateInput = p.Results.EndDate;

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

raw = load(matFile);
if ~isfield(raw, 'datos_96')
    error('run_adaptivity_figures_datos96_unified:MissingTable', ...
        'El archivo %s no contiene la tabla datos_96.', matFile);
end

if ~istable(raw.datos_96)
    error('run_adaptivity_figures_datos96_unified:InvalidTable', ...
        'datos_96 no es una tabla válida.');
end

tbl = raw.datos_96;
if ~ismember('timestamp', tbl.Properties.VariableNames)
    error('run_adaptivity_figures_datos96_unified:MissingTimestamp', ...
        'La tabla datos_96 no contiene la columna timestamp.');
end

% Parse timestamp to datetime and sort

t = datos96_parse_timestamp(tbl.timestamp);
[tsorted, order] = sort(t);
tbl = tbl(order, :);

% Remove duplicate timestamps, keep last occurrence
[tu, idxLast] = unique(tsorted, 'last');
tbl = tbl(idxLast, :);
t = tu;

% Filter by date range (inclusive). Default: 2024-01-12 to 2024-02-05.
if isempty(startDateInput)
    startDateInput = '2024-01-12 00:00:00';
end
if isempty(endDateInput)
    endDateInput = '2024-02-05 23:59:59';
end
startDate = datos96_parse_date_input(startDateInput, t);
endDate = datos96_parse_date_input(endDateInput, t);

if ~isempty(startDate) || ~isempty(endDate)
    preCount = height(tbl);
    if isempty(startDate)
        startDate = min(t);
    end
    if isempty(endDate)
        endDate = max(t);
    end
    mask = (t >= startDate) & (t <= endDate);
    tbl = tbl(mask, :);
    t = t(mask);
    postCount = height(tbl);
    if postCount == 0
        error('run_adaptivity_figures_datos96_unified:EmptyDateRange', ...
            ['El filtro de fechas dejó la tabla vacía. ' ...
             'Filas antes: %d, filas después: %d.'], preCount, postCount);
    end
end

% Extract variables (degrade gracefully if missing)
A11 = datos96_get_table_var(tbl, 'A11');
A12 = datos96_get_table_var(tbl, 'A12');
A21 = datos96_get_table_var(tbl, 'A21');
A22 = datos96_get_table_var(tbl, 'A22');
B1 = datos96_get_table_var(tbl, 'B1');
B2 = datos96_get_table_var(tbl, 'B2');

if ~ismember('uk_1', tbl.Properties.VariableNames)
    error('run_adaptivity_figures_datos96_unified:MissingUk1', ...
        'La tabla datos_96 no contiene la columna obligatoria uk_1.');
end
u_applied = datos96_get_table_var(tbl, 'uk_1');
u_mpc = datos96_get_table_var(tbl, 'uk_mpc');

% Build theta vector with NaN placeholders for missing columns
n = height(tbl);
theta = nan(n, 6);
theta(:, 1) = A11;
theta(:, 2) = A21;
theta(:, 3) = A12;
theta(:, 4) = A22;
theta(:, 5) = B1;
theta(:, 6) = B2;

% Delta theta and norm

dtheta = [nan(1, size(theta, 2)); diff(theta, 1, 1)];
dtheta_norm = sqrt(sum(dtheta.^2, 2, 'omitnan'));
idxAllNaN = all(isnan(dtheta), 2);
dtheta_norm(idxAllNaN) = NaN;

% Stats
validDtheta = dtheta_norm(~isnan(dtheta_norm));
if isempty(validDtheta)
    stats.dtheta_mean = NaN;
    stats.dtheta_p95 = NaN;
    stats.dtheta_p99 = NaN;
else
    stats.dtheta_mean = mean(validDtheta);
    stats.dtheta_p95 = prctile(validDtheta, 95);
    stats.dtheta_p99 = prctile(validDtheta, 99);
end

u_on = (u_applied > 0.5);
u_on(isnan(u_applied)) = false;
stats.percent_irrigation_on = mean(u_on) * 100;
% "irrigation_events" cuenta flancos ascendentes (0 -> 1) en u_on.
stats.irrigation_events = sum(diff(u_on) == 1);

% Helper for semilog plotting of ||Delta theta|| (avoid zeros/negatives)
validPos = validDtheta(validDtheta > 0);
if isempty(validPos)
    dthetaFloor = 1e-10;
else
    dthetaFloor = max(1e-10, min(validPos) / 10);
end
dtheta_plot = dtheta_norm;
dtheta_plot(~isnan(dtheta_plot) & dtheta_plot <= 0) = dthetaFloor;

% Prepare output struct
S = struct();
S.t = t;
S.u = u_applied;
S.u_mpc = u_mpc;
S.A11 = A11;
S.A12 = A12;
S.A21 = A21;
S.A22 = A22;
S.B1 = B1;
S.B2 = B2;
S.theta = theta;
S.dtheta = dtheta;
S.dtheta_norm = dtheta_norm;

% Figure 1: A parameters
figA = figure('Color', 'w');
hold on;
if ~all(isnan(A11))
    plot(t, A11, 'Color', [0.1 0.1 0.5], 'LineWidth', 1.2, ...
        'DisplayName', 'A_{11}');
end
if ~all(isnan(A22))
    plot(t, A22, 'Color', [0.5 0.1 0.1], 'LineWidth', 1.2, ...
        'DisplayName', 'A_{22}');
end
if ~all(isnan(A12))
    plot(t, A12, 'Color', [0.2 0.5 0.2], 'LineWidth', 1.1, ...
        'DisplayName', 'A_{12}');
end
if ~all(isnan(A21))
    plot(t, A21, 'Color', [0.4 0.2 0.6], 'LineWidth', 1.1, ...
        'DisplayName', 'A_{21}');
end
hold off;

xlabel('Time');
ylabel('A-matrix parameters');
title('Temporal evolution of A-matrix parameters (RLS)');
legend('Location', 'eastoutside');
grid on;

if all(isnan(A11)) && all(isnan(A22)) && all(isnan(A12)) && all(isnan(A21))
    warning('No hay parámetros A disponibles para graficar.');
end

figAPathPng = fullfile(saveDir, 'fig01_parametros_A.png');
figAPathPdf = fullfile(saveDir, 'fig01_parametros_A.pdf');
datos96_export_figure(figA, figAPathPng, figAPathPdf);

% Figure 2: B parameters
figB = figure('Color', 'w');
hold on;
if ~all(isnan(B1))
    plot(t, B1, 'Color', [0 0.3 0.6], 'LineWidth', 1.2, ...
        'DisplayName', 'B_{1}');
end
if ~all(isnan(B2))
    plot(t, B2, 'Color', [0.6 0.3 0], 'LineWidth', 1.2, ...
        'DisplayName', 'B_{2}');
end
hold off;

xlabel('Time');
ylabel('B-matrix parameters');
title('Temporal evolution of B-matrix parameters (RLS)');
grid on;

if shadeIrrigation && ~isempty(u_applied)
    datos96_add_irrigation_shading(figB, t, u_applied);
end

legend('Location', 'eastoutside');

if all(isnan(B1)) && all(isnan(B2))
    warning('No hay parámetros B disponibles para graficar.');
end

figBPathPng = fullfile(saveDir, 'fig02_parametros_B.png');
figBPathPdf = fullfile(saveDir, 'fig02_parametros_B.pdf');
datos96_export_figure(figB, figBPathPng, figBPathPdf);

% Figure 3: ||Delta theta||_2
figD = figure('Color', 'w');
if useSemilogDtheta
    plot(t, dtheta_plot, 'Color', [0.2 0.2 0.7], 'DisplayName', '||\Delta\theta||_2');
    set(gca, 'YScale', 'log');
else
    plot(t, dtheta_norm, 'Color', [0.2 0.2 0.7], 'DisplayName', '||\Delta\theta||_2');
end

hold on;
try
    tt = timetable(t, dtheta_norm);
    switch dailyAgg
        case 'mean'
            ttDay = retime(tt, 'daily', 'mean');
            aggName = 'Daily mean';
        otherwise
            ttDay = retime(tt, 'daily', 'median');
            aggName = 'Daily median';
    end
    yAgg = ttDay.dtheta_norm;
    if useSemilogDtheta
        yAgg(~isnan(yAgg) & yAgg <= 0) = dthetaFloor;
    end
    plot(ttDay.t, yAgg, '-o', 'Color', [0.8 0.2 0.2], ...
        'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', aggName);
catch
    warning('No se pudo calcular la agregación diaria de ||Delta theta||.');
end
hold off;

xlabel('Time');
if useSemilogDtheta
    ylabel('||\Delta\theta||_2 (log scale)');
else
    ylabel('||\Delta\theta||_2');
end
title('Temporal evolution of adaptation magnitude');
legend('Location', 'northwest');
grid on;

if shadeIrrigation && ~isempty(u_applied)
    datos96_add_irrigation_shading(figD, t, u_applied);
end

figDPathPng = fullfile(saveDir, 'fig03_dtheta_norm.png');
figDPathPdf = fullfile(saveDir, 'fig03_dtheta_norm.pdf');
datos96_export_figure(figD, figDPathPng, figDPathPdf);


% Figure 4: Unified adaptivity overview (A, B, and ||Delta theta||_2)
figU = figure('Color', 'w');
tl = tiledlayout(figU, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Field demonstration of adaptivity (RLS): A(k), B(k), and ||\Delta\theta(k)||_2', 'FontSize', 16);

% --- Panel (a): A-matrix parameters ---
ax1 = nexttile(tl, 1);
hold(ax1, 'on');
if ~all(isnan(A11)), plot(ax1, t, A11, 'Color', [0.1 0.1 0.5], 'LineWidth', 1.2, 'DisplayName', 'A_{11}'); end
if ~all(isnan(A22)), plot(ax1, t, A22, 'Color', [0.5 0.1 0.1], 'LineWidth', 1.2, 'DisplayName', 'A_{22}'); end
if ~all(isnan(A12)), plot(ax1, t, A12, 'Color', [0.2 0.5 0.2], 'LineWidth', 1.1, 'DisplayName', 'A_{12}'); end
if ~all(isnan(A21)), plot(ax1, t, A21, 'Color', [0.4 0.2 0.6], 'LineWidth', 1.1, 'DisplayName', 'A_{21}'); end
hold(ax1, 'off');
ylabel(ax1, 'A-matrix entries');
grid(ax1, 'on');
legend(ax1, 'Location', 'eastoutside');
ax1.XTickLabel = []; % keep x-axis labels only on the bottom panel

if shadeIrrigation && ~isempty(u_applied)
    axes(ax1); %#ok<LAXES>
    datos96_add_irrigation_shading(figU, t, u_applied);
end

% --- Panel (b): B-matrix parameters (with irrigation shading) ---
ax2 = nexttile(tl, 2);
hold(ax2, 'on');
if ~all(isnan(B1)), plot(ax2, t, B1, 'Color', [0 0.3 0.6], 'LineWidth', 1.2, 'DisplayName', 'B_{1}'); end
if ~all(isnan(B2)), plot(ax2, t, B2, 'Color', [0.6 0.3 0], 'LineWidth', 1.2, 'DisplayName', 'B_{2}'); end
hold(ax2, 'off');
ylabel(ax2, 'B-matrix entries');
grid(ax2, 'on');
if shadeIrrigation && ~isempty(u_applied)
    % Shade irrigation ON periods based on logged applied valve state uk_1
    axes(ax2); %#ok<LAXES>
    datos96_add_irrigation_shading(figU, t, u_applied);
end
legend(ax2, 'Location', 'eastoutside');
ax2.XTickLabel = [];

% --- Panel (c): ||Delta theta||_2 (with daily aggregation and irrigation shading) ---
ax3 = nexttile(tl, 3);
if useSemilogDtheta
    plot(ax3, t, dtheta_plot, 'Color', [0.2 0.2 0.7], 'DisplayName', '||\Delta\theta||_2');
    set(ax3, 'YScale', 'log');
else
    plot(ax3, t, dtheta_norm, 'Color', [0.2 0.2 0.7], 'DisplayName', '||\Delta\theta||_2');
end
hold(ax3, 'on');
try
    tt = timetable(t, dtheta_norm);
    switch dailyAgg
        case 'mean'
            ttDay = retime(tt, 'daily', 'mean');
            aggName = 'Daily mean';
        otherwise
            ttDay = retime(tt, 'daily', 'median');
            aggName = 'Daily median';
    end
    yAgg = ttDay.dtheta_norm;
    if useSemilogDtheta
        yAgg(~isnan(yAgg) & yAgg <= 0) = dthetaFloor;
    end
    plot(ax3, ttDay.t, yAgg, '-o', 'Color', [0.8 0.2 0.2], ...
        'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', aggName);
catch
    warning('No se pudo calcular la agregación diaria de ||Delta theta||.');
end
hold(ax3, 'off');
xlabel(ax3, 'Time');
if useSemilogDtheta
    ylabel(ax3, '||\Delta\theta||_2 (log scale)');
else
    ylabel(ax3, '||\Delta\theta||_2');
end
grid(ax3, 'on');
if shadeIrrigation && ~isempty(u_applied)
    axes(ax3); %#ok<LAXES>
    datos96_add_irrigation_shading(figU, t, u_applied);
end
legend(ax3, 'Location', 'northwest');

% Link x-axes across panels
linkaxes([ax1, ax2, ax3], 'x');

% Export unified figure
figUPathPng = fullfile(saveDir, 'fig04_adaptivity_unified.png');
figUPathPdf = fullfile(saveDir, 'fig04_adaptivity_unified.pdf');
datos96_export_figure(figU, figUPathPng, figUPathPdf);

end
% [S, stats] = run_adaptivity_figures_datos96_unified('datos_96_act.mat', ...
%     'SaveDir','figs', ...
%     'DailyAgg','median', ...
%     'StartDate','2024-01-12 00:00:00', ...
%     'EndDate',  '2024-02-05 23:59:59');

% [S, stats] = run_adaptivity_figures_datos96_unified('datos_96_act.mat', ...
%     'SaveDir','figs', ...
%     'ShadeIrrigation', false);

