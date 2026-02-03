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
fecha_inicio = datetime('2024-01-16 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
fecha_fin = datetime('2024-02-01 23:59:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Filtrar los datos según la opción seleccionada
if ~graficar_todos
    data = data(data.timestamp96 >= fecha_inicio & data.timestamp96 <= fecha_fin, :);
end

%% gráfica de humedades
figure(1);
% subplot(2,1,1);
ylim([0.5 0.9]);
hold on;

% Graficar las series de datos con formatos específicos
plot(data.timestamp96, data.x1k, 'r.', 'MarkerSize', 5); % x1
plot(data.timestamp96, data.x2k, 'b.', 'MarkerSize', 5); % x2
plot(data.timestamp96, data.x1min, 'r--', 'LineWidth', 1); % x1min
plot(data.timestamp96, data.x1max, 'r--', 'LineWidth', 1); % x1max
plot(data.timestamp96, data.x2min, 'b--', 'LineWidth', 1); % x2min
plot(data.timestamp96, data.x2max, 'b--', 'LineWidth', 1); % x2max

% Configuración de etiquetas y leyenda
legend({'x1', 'x2', 'x1min', 'x1max', 'x2min', 'x2max'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Humidity');
%title('Series Temporales de Humedad');
hold off;

%% gráfica de parámetros
% Configuración de la figura y subgráficos
figure(2);
% subplot(2,1,1);
%ylim([0.05 0.1]);
hold on;

% Graficar las series de datos con formatos específicos
% plot(data.timestamp96, data.A11, 'r--', 'MarkerSize', 2); % x1
% plot(data.timestamp96, data.A12, 'g--', 'MarkerSize', 2); % x2
plot(data.timestamp96, data.A21, 'b.', 'LineWidth', 2); % x1min
%plot(data.timestamp96, data.A22, 'b.', 'LineWidth', 2); % x1max
% plot(data.timestamp96, data.B1, 'k--', 'LineWidth', 2); % x2min
% plot(data.timestamp96, data.B2, 'c--', 'LineWidth', 2); % x2max

% Configuración de etiquetas y leyenda
% legend({'A11', 'A12', 'A21', 'A22', 'B1', 'B2'}, 'Location', 'best', 'Orientation', 'horizontal');
legend({'A21'}, 'Location', 'best', 'Orientation', 'horizontal');
xlabel('Time');
ylabel('Paremeter');
%title('Series Temporales de Humedad');
hold off;

% Subgráfico para la serie uk_mpc96
% subplot(2,1,2);
% 
% stairs(data.timestamp96, data.uk_mpc96, 'k-', 'LineWidth', 1.5); % uk_mpc96
% ylim([0 1.1]);
% xlabel('Tiempo');
% ylabel('Control MPC');
% title('Control MPC a lo largo del tiempo');

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