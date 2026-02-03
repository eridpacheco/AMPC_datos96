clear all
clc
close all

% Cargar los datos y convertir la columna de tiempo a formato datetime
load('datos_96_act.mat');
data = datos_96;
data.timestamp96 = datetime(data.timestamp, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Parámetro para graficar todos los datos o solo un rango de fechas
graficar_todos = false; % Cambiar a true para graficar todos los datos

% Rango de fechas en caso de seleccionar un subconjunto
% fecha_inicio = datetime('2024-01-13 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
% fecha_fin = datetime('2024-02-05 23:59:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
fecha_inicio = datetime('2024-01-15 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
fecha_fin    = datetime('2024-01-15 23:59:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Filtrar los datos según la opción seleccionada
if ~graficar_todos
    data = data(data.timestamp96 >= fecha_inicio & data.timestamp96 <= fecha_fin, :);
end

%% gráfica de humedades
figure(1);
subplot(2,1,1);
ylim([0.5 0.9]);
hold on;

% Graficar las series de datos con formatos específicos
%plot(data.timestamp96, data.x1k, 'r.', 'MarkerSize', 5); % x1
plot(data.timestamp96, data.x2k, 'b.', 'MarkerSize', 5); % x2
%plot(data.timestamp96, data.x1min, 'r--', 'LineWidth', 1); % x1min
%plot(data.timestamp96, data.x1max, 'r--', 'LineWidth', 1); % x1max
plot(data.timestamp96, data.x2min, 'b--', 'LineWidth', 1); % x2min
plot(data.timestamp96, data.x2max, 'b--', 'LineWidth', 1); % x2max

% Configuración de etiquetas y leyenda
legend({'x2', 'x2min', 'x2max'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Measured Humidity of Layer 2');
%title('Series Temporales de Humedad');
hold off;

subplot(2,1,2);
ylim([0.5 0.9]);
hold on;

% Graficar las series de datos con formatos específicos
%plot(data.timestamp96, data(1,:).x1_pred, 'r.', 'MarkerSize', 5); % x1
plot(data.timestamp96, data(1,:).x2_pred', 'b.', 'MarkerSize', 5); % x2
%plot(data.timestamp96, data.x1min, 'r--', 'LineWidth', 1); % x1min
%plot(data.timestamp96, data.x1max, 'r--', 'LineWidth', 1); % x1max
plot(data.timestamp96, data.x2min, 'b--', 'LineWidth', 1); % x2min
plot(data.timestamp96, data.x2max, 'b--', 'LineWidth', 1); % x2max

% Configuración de etiquetas y leyenda
legend({'x2_{pred}', 'x2min', 'x2max'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Predicted Humidity of Layer 2');
%title('Series Temporales de Humedad');
hold off;
a1 = data.x2k;
a2 = data(1,:).x2_pred';
error = rmse(a1,a2);

%% gráfica de parámetros
% Configuración de la figura y subgráficos
figure(2);
subplot(3,1,1);
%ylim([0.05 0.1]);
hold on;
plot(data.timestamp96, data.A21, 'k.', 'LineWidth', 2); % x1min
legend({'A21'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Paremeter');
hold off;

subplot(3,1,2);
%ylim([0.05 0.1]);
hold on;
plot(data.timestamp96, data.A22, 'k.', 'LineWidth', 2); % x1min
legend({'A22'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Paremeter');
hold off;

subplot(3,1,3);
%ylim([0.05 0.1]);
hold on;
plot(data.timestamp96, data.B2, 'k.', 'LineWidth', 2); % x1min
legend({'B2'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Paremeter');
hold off;

% Cálculos adicionales
cantidad_unos = sum(data.uk_mpc);
nro_dias = days(fecha_fin - fecha_inicio);
prom_diario = cantidad_unos / nro_dias;
vol_total = (cantidad_unos * 15 / 60) * 1.6 * 100 * 4;
vol_diario = vol_total / nro_dias;

% Mostrar los resultados en la consola
fprintf('Cantidad de unos en uk_mpc96: %d\n', cantidad_unos);
fprintf('Número de días: %d\n', nro_dias);
fprintf('Promedio diario de unos: %.2f\n', prom_diario);
fprintf('Volumen total: %.2f\n', vol_total);
fprintf('Volumen diario: %.2f\n', vol_diario);
disp("Error = " + error)

