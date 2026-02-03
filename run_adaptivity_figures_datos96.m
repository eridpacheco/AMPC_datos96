function [S, stats] = run_adaptivity_figures_datos96(matFile, varargin)
%RUN_ADAPTIVITY_FIGURES_DATOS96 Generate adaptivity figures for datos_96.


%[S, stats] = run_adaptivity_figures_datos96('datos_96_act.mat', ...
%    'SaveDir', 'figs', ...
%    'WindowSamples', 96, ...
%    'DailyAgg', 'median');

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
% Smoothing uses a centered moving median (default window = 96 samples,
% corresponding to 1 day at 15-min sampling):
%   smooth_x(k) = movmedian(x, windowSamples, 'omitnan').
%
% Usage:
%   [S, stats] = run_adaptivity_figures_datos96('datos_96_act.mat', ...
%       'SaveDir', 'figs', 'WindowSamples', 96, 'DailyAgg', 'median');

p = inputParser;
p.addRequired('matFile', @(x) ischar(x) || isstring(x));
p.addParameter('SaveDir', 'figs', @(x) ischar(x) || isstring(x));
p.addParameter('WindowSamples', 96, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('DailyAgg', 'median', @(x) ischar(x) || isstring(x));
p.parse(matFile, varargin{:});

saveDir = char(p.Results.SaveDir);
windowSamples = p.Results.WindowSamples;
dailyAgg = lower(char(p.Results.DailyAgg));

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

raw = load(matFile);
if ~isfield(raw, 'datos_96')
    error('run_adaptivity_figures_datos96:MissingTable', ...
        'El archivo %s no contiene la tabla datos_96.', matFile);
end

if ~istable(raw.datos_96)
    error('run_adaptivity_figures_datos96:InvalidTable', ...
        'datos_96 no es una tabla válida.');
end

tbl = raw.datos_96;
if ~ismember('timestamp', tbl.Properties.VariableNames)
    error('run_adaptivity_figures_datos96:MissingTimestamp', ...
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

% Extract variables (degrade gracefully if missing)
A11 = datos96_get_table_var(tbl, 'A11');
A12 = datos96_get_table_var(tbl, 'A12');
A21 = datos96_get_table_var(tbl, 'A21');
A22 = datos96_get_table_var(tbl, 'A22');
B1 = datos96_get_table_var(tbl, 'B1');
B2 = datos96_get_table_var(tbl, 'B2');

u = datos96_get_table_var(tbl, 'uk_1');

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

if isempty(u)
    stats.percent_irrigation_on = NaN;
    stats.irrigation_events = NaN;
else
    u_on = (u == 1);
    stats.percent_irrigation_on = mean(u_on, 'omitnan') * 100;
    stats.irrigation_events = sum(diff(u_on) == 1, 'omitnan');
end

% Prepare output struct
S = struct();
S.t = t;
S.u = u;
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
plotA = {};
if ~all(isnan(A11))
    plot(t, A11, 'Color', [0.6 0.6 0.9], 'DisplayName', 'A11 (crudo)');
    plot(t, movmedian(A11, windowSamples, 'omitnan'), 'Color', [0.1 0.1 0.5], ...
        'LineWidth', 1.6, 'DisplayName', 'A11 (suavizado)');
    plotA{end+1} = 'A11';
end
if ~all(isnan(A22))
    plot(t, A22, 'Color', [0.9 0.6 0.6], 'DisplayName', 'A22 (crudo)');
    plot(t, movmedian(A22, windowSamples, 'omitnan'), 'Color', [0.5 0.1 0.1], ...
        'LineWidth', 1.6, 'DisplayName', 'A22 (suavizado)');
    plotA{end+1} = 'A22';
end
if ~all(isnan(A12))
    plot(t, A12, 'Color', [0.7 0.8 0.7], 'DisplayName', 'A12 (crudo)');
    plot(t, movmedian(A12, windowSamples, 'omitnan'), 'Color', [0.2 0.5 0.2], ...
        'LineWidth', 1.4, 'DisplayName', 'A12 (suavizado)');
    plotA{end+1} = 'A12';
end
if ~all(isnan(A21))
    plot(t, A21, 'Color', [0.8 0.7 0.9], 'DisplayName', 'A21 (crudo)');
    plot(t, movmedian(A21, windowSamples, 'omitnan'), 'Color', [0.4 0.2 0.6], ...
        'LineWidth', 1.4, 'DisplayName', 'A21 (suavizado)');
    plotA{end+1} = 'A21';
end
hold off;

xlabel('Tiempo');
ylabel('Parámetros A');
title('Evolución temporal de parámetros A (RLS)');
legend('Location', 'best');
grid on;

if isempty(plotA)
    warning('No hay parámetros A disponibles para graficar.');
end

figAPathPng = fullfile(saveDir, 'fig01_parametros_A.png');
figAPathPdf = fullfile(saveDir, 'fig01_parametros_A.pdf');
datos96_export_figure(figA, figAPathPng, figAPathPdf);

% Figure 2: B parameters
figB = figure('Color', 'w');
hold on;
plotB = {};
if ~all(isnan(B1))
    plot(t, B1, 'Color', [0.5 0.7 0.9], 'DisplayName', 'B1 (crudo)');
    plot(t, movmedian(B1, windowSamples, 'omitnan'), 'Color', [0 0.3 0.6], ...
        'LineWidth', 1.6, 'DisplayName', 'B1 (suavizado)');
    plotB{end+1} = 'B1';
end
if ~all(isnan(B2))
    plot(t, B2, 'Color', [0.9 0.7 0.5], 'DisplayName', 'B2 (crudo)');
    plot(t, movmedian(B2, windowSamples, 'omitnan'), 'Color', [0.6 0.3 0], ...
        'LineWidth', 1.6, 'DisplayName', 'B2 (suavizado)');
    plotB{end+1} = 'B2';
end
hold off;

xlabel('Tiempo');
ylabel('Parámetros B');
title('Evolución temporal de parámetros B (RLS)');
grid on;

if ~isempty(u)
    datos96_add_irrigation_shading(figB, t, u);
end

legend('Location', 'best');

if isempty(plotB)
    warning('No hay parámetros B disponibles para graficar.');
end

figBPathPng = fullfile(saveDir, 'fig02_parametros_B.png');
figBPathPdf = fullfile(saveDir, 'fig02_parametros_B.pdf');
datos96_export_figure(figB, figBPathPng, figBPathPdf);

% Figure 3: ||Delta theta||_2
figD = figure('Color', 'w');
plot(t, dtheta_norm, 'Color', [0.2 0.2 0.7], 'DisplayName', '||\Delta\theta||_2');

hold on;
try
    tt = timetable(t, dtheta_norm);
    switch dailyAgg
        case 'mean'
            ttDay = retime(tt, 'daily', 'mean');
            aggName = 'Media diaria';
        otherwise
            ttDay = retime(tt, 'daily', 'median');
            aggName = 'Mediana diaria';
    end
    plot(ttDay.t, ttDay.dtheta_norm, '-o', 'Color', [0.8 0.2 0.2], ...
        'LineWidth', 1.2, 'MarkerSize', 4, 'DisplayName', aggName);
catch
    warning('No se pudo calcular la agregación diaria de ||Delta theta||.');
end
hold off;

xlabel('Tiempo');
ylabel('||\Delta\theta||_2');
title('Evolución de la magnitud de adaptación');
legend('Location', 'best');
grid on;

figDPathPng = fullfile(saveDir, 'fig03_dtheta_norm.png');
figDPathPdf = fullfile(saveDir, 'fig03_dtheta_norm.pdf');
datos96_export_figure(figD, figDPathPng, figDPathPdf);

end