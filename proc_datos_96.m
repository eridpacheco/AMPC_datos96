% Define el directorio donde están los archivos .mat
folder_path = 'C:\Users\eridp\OneDrive - Universidad Loyola Andalucía\Documentos\Tesis\Papers\A enviar\Revista AMPC\Datos\96';

% Obtén una lista de todos los archivos .mat en el directorio
files = dir(fullfile(folder_path, '*.mat'));

% Variables para almacenar los datos concatenados
timestamp = [];
controller = [];
sname = [];
error = [];
x1min = [];
x2min = [];
x1max = [];
x2max = [];
x1k = [];
x2k = [];
uk_1 = [];
uk_mpc = [];
x1_pred = [];
x2_pred = [];
k_mpc = [];
u_pred = [];
z_cost = [];
A11 = [];
A12 = [];
A21 = [];
A22 = [];
B1 = [];
B2 = [];
c1 = 0;
ind_ij = [];
% Itera a través de cada archivo
for i = 1:length(files)
    % Carga el archivo .mat actual
    file_path = fullfile(folder_path, files(i).name);
    load(file_path);
    disp(file_path)
    for j = 1:length(MPCvar1)
        try
            if MPCvar1(j).Controller == 5
                x1_pred = [x1_pred; MPCvar1(j).mpc.x1_sal];
                x2_pred = [x2_pred; MPCvar1(j).mpc.x2_sal];
                k_mpc = [k_mpc; MPCvar1(j).mpc.k];
                u_pred = [u_pred; MPCvar1(j).mpc.u_sal];
                z_cost = [z_cost; MPCvar1(j).mpc.z_sal];
                A11 = [A11; MPCvar1(j).model.A11];
                A12 = [A12; MPCvar1(j).model.A12];
                A21 = [A21; MPCvar1(j).model.A21];
                A22 = [A22; MPCvar1(j).model.A22];
                B1 = [B1; MPCvar1(j).model.B1];
                B2 = [B2; MPCvar1(j).model.B2];
                x1min = [x1min; MPCvar1(j).mpc.x1min];
                x2min = [x2min; MPCvar1(j).mpc.x2min];
                x1max = [x1max; MPCvar1(j).mpc.x1max];
                x2max = [x2max; MPCvar1(j).mpc.x2max];
                timestamp = [timestamp; MPCvar1(j).timek];
                controller = [controller; MPCvar1(j).Controller];
                sname = [sname; MPCvar1(j).sname];
                x1k = [x1k; MPCvar1(j).x1k_e];
                x2k = [x2k; MPCvar1(j).x2k_e];
                uk_1 = [uk_1; MPCvar1(j).uk_1];
                uk_mpc = [uk_mpc; MPCvar1(j).mpc.uk];
            end
        catch except
            c1 = c1 + 1;
            a = [file_path string(j)];
            ind_ij = [ind_ij; a];
        end
    end
end
% Crea una tabla con los datos
datos_96 = table(timestamp, controller, k_mpc, sname, x1k, x2k, uk_1, uk_mpc, x1_pred, x2_pred,...
                u_pred, z_cost, A11, A12, A21, A22, B1, B2, x1min, x2min, x1max,x2max);
% Guarda los datos concatenados en un nuevo archivo .mat
save(fullfile(folder_path, 'datos_96_act.mat'), 'datos_96');
