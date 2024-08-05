file_name = "COMSOL\\RibbonSynapse";
component_name = "RibbonSynapse";

interpolation_unit = 's';

cell_vesicle_distance = 60;

import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

% clear model;

if ~exist("model", "var")
    try
       model = mphload(file_name);
    catch E
       fprintf("Please create the RibbonSynapse file %s first\n", file_name);
       return;
    end
end

try
   model.component(component_name);
catch
   fprintf("Please define the model component %s first\n", component_name);
   return;
end

labels = {"r", "m", "i1", "i2", "i3", "i4", "i5"};

for j = 1:length(labels)
    rng(j + 21, "twister");
    
    label = labels{j};
    rs = getSynapse(label);

    opts = SynapseOptions;
    num_time_points = 1e4;
    dt = 1e-4;
    % Vt = -0.07;
    to_csv = false;
    
    create_geometry(model, rs, component_name, cell_vesicle_distance);
    create_channel_selections(model, rs, component_name);
    create_vesicle_selections(model, rs, component_name);
    create_integrations(model, rs, component_name);

    top_side = setdiff(model.selection("box3").entities, model.selection("box1").entities);
    tds = model.physics("tds");
    tds.selection.set(1)
    tds.feature("CalciumConcentration").selection.set(top_side);
    
    voltage_steps = linspace(-.060, -.00, 40);
    
    global_eval = model.result.numerical("gev1");
    
    disp("Voltage list"); disp(voltage_steps);
    
    for i = 1:40

        Vt = voltage_steps(i);
        Vt_STR = num2str(abs(Vt));

        fprintf("Running Voltage Settings on Vt %s\n", Vt_STR(3:end));

        model.param("par9").set("Vm", sprintf("%d [V]", Vt));

        channel_state_matrix = GenerateChannelUpdateMatrix(opts, 1e4, 1e-4, Vt);
        create_flux("block", component_name, rs, model, channel_state_matrix)

        try
            model.study("std2").run;
        catch E
            disp(E)
        end

        global_eval.setResult();

        results = mphtable(model, "tbl1");
        tbl_data = results.data;

        [status, msg, msgID] = mkdir(sprintf("Data-%s", label));
        writematrix(tbl_data, sprintf("Data-%s\\Vt_%s.csv", label, Vt_STR(3:end)));
    end

end


%% Update Geometry
function create_geometry(model, rs, component_name, cell_vesicle_distance)
fig = model.component(component_name).geom('Geometry');

% Clear existing geometry features without removing the entire geometry
tags = fig.feature.tags;
for i = 1:length(tags)
    if tags(i).toCharArray' == "fin"
        continue
    end
    fig.feature.remove(tags(i));
end

% Redefine the geometry features
max_x = max([rs.xv(:, 1); rs.xc(:, 1)]) + rs.vesicle_radius + cell_vesicle_distance;
min_x = min([rs.xv(:, 1); rs.xc(:, 1)]) - rs.vesicle_radius - cell_vesicle_distance;
max_y = max([rs.xv(:, 2); rs.xc(:, 2)]) + rs.vesicle_radius + cell_vesicle_distance;
min_y = min([rs.xv(:, 2); rs.xc(:, 2)]) - rs.vesicle_radius - cell_vesicle_distance;
max_z = max([rs.xv(:, 3); rs.xc(:, 3)]) + rs.vesicle_radius + cell_vesicle_distance;

%% Define Cell
blockName = 'CellBlock';
fig.feature.create(blockName, 'Block');
fig.feature(blockName).set('size', [max_x - min_x, max_y - min_y, max_z]);
fig.feature(blockName).set('pos', [min_x, min_y, 0]);

%% Define Vesicles
for i = 1:size(rs.xv, 1)
    sphereName = ['Vesicle_' num2str(i)];
    fig.feature.create(sphereName, 'Sphere');
    fig.feature(sphereName).set('r', rs.vesicle_radius);
    fig.feature(sphereName).set('pos', [rs.xv(i, 1), rs.xv(i, 2), rs.xv(i, 3)]);
    fig.feature(sphereName).label(sprintf("Vesicle %i", i));

end

%% Define Channels
wpName = 'WorkPlaneChannels';
fig.feature.create(wpName, 'WorkPlane');
fig.feature(wpName).set('quickplane', 'xy');
fig.feature(wpName).set('quickz', rs.xc(1, 3)); % Assuming all channels are in the same z-plane
for i = 1:size(rs.xc, 1)
    circleName = ['Channel_' num2str(i)];
    fig.feature(wpName).geom.feature.create(circleName, 'Circle');
    fig.feature(wpName).geom.feature(circleName).set('r', rs.channel_radius);
    fig.feature(wpName).geom.feature(circleName).set('pos', [rs.xc(i, 1), rs.xc(i, 2)]);
    fig.feature(wpName).geom.feature(circleName).label(sprintf("Channel %i", i));
end

fig.run;

try
    model.component(component_name).mesh.remove("Normal_Mesh");
end

mesh = model.component(component_name).mesh.create("Normal_Mesh");
mesh.label("Normal");

end


%% Create Channel Selections
function create_channel_selections(model, rs, component_name)
channel_indices = size(rs.xc, 1);

for i = 1:size(rs.xc, 1)
    try
        ch_sel = model.component(component_name).selection.create(sprintf("Channel_Sel_%i", i), "Ball");
    catch
        ch_sel = model.component(component_name).selection(sprintf("Channel_Sel_%i", i));
    end
    ch_sel.set("posx", rs.xc(i, 1)).set("posy", rs.xc(i, 2)).set("posz", rs.xc(i, 3));
    ch_sel.set("r", 0);
    ch_sel.set("entitydim", 2);
    ch_sel.label(sprintf("Channel %i", i))


    channel_indices(i) = ch_sel.entities;
end

try
    ch_sel = model.component(component_name).selection.create("Channels_Sel", "Explicit");
catch
    ch_sel = model.component(component_name).selection("Channels_Sel");
end
ch_sel.geom(2);
ch_sel.set(channel_indices);
ch_sel.label("Channels");
end

%% Create Vesicle Selections
function create_vesicle_selections(model, rs, component_name)
for i = 1:size(rs.xv, 1)
    try
        ves_sel = model.component(component_name).selection.create(sprintf("Vesicle_Sel_%i", i), "Ball");
    catch
        ves_sel = model.component(component_name).selection(sprintf("Vesicle_Sel_%i", i));
    end
    ves_sel.set("posx", rs.xv(i, 1)).set("posy", rs.xv(i, 2)).set("posz", rs.xv(i, 3));
    ves_sel.set("r", 20);
    ves_sel.set("entitydim", 2);
    ves_sel.label(sprintf("Vesicle %i", i));
end
end

%% Create Integrations
function create_integrations(model, rs, component_name)
for i = 1:size(rs.xv, 1)
    ves_sel = model.component(component_name).selection(sprintf("Vesicle_Sel_%i", i));

    cpl = model.cpl(sprintf("intop%i", i));
    % cpl = model.cpl.create(["intop" i], "Integration");
    % cpl.geom("geom", 2);
    cpl.selection.set(ves_sel.entities);

    model.component(component_name).variable("Average Concentrations").set(sprintf("C_Avg_%i", i), sprintf("intop%i(c)/intop%i(1)", i, i)).label("Average Concentrations");
end
end

%% Create Fluxes
function create_flux(setting, component_name, rs, model, channel_state_matrix)
arguments
    setting string;
    component_name string;
    rs;
    model;
    channel_state_matrix (:, :) double = [];
end

channel_indices = size(rs.xc, 1);
tds = model.physics("tds");

% tags = model.component(component_name').selection.tags;
% for i = 1:length(tags)
%     if ~startsWith(tags(i).toCharArray', "Channel_Sel_")
%         continue
%     end
%     model.component(component_name').selection.remove(tags(i))
% end
% 
% for i = 1:size(rs.xc, 1)
%     ch_sel = model.component(component_name').selection.create(sprintf("Channel_Sel_%i", i), "Ball");
%     ch_sel.set("posx", rs.xc(i, 1)).set("posy", rs.xc(i, 2)).set("posz", rs.xc(i, 3));
%     ch_sel.set("r", 0);
%     ch_sel.set("entitydim", 2);
%     channel_indices(i) = ch_sel.entities;
% end

if setting == "complex" | setting == "block"
    % tags = tds.feature.tags;
    % for i = 1:length(tags)
    %     if ~startsWith(tags(i).toCharArray', "Flux_Ch_")
    %         continue
    %     end
    %     tds.feature.remove(tags(i))
    % end

    for i = 1:size(rs.xc, 1)
        ch_sel = model.component(component_name').selection(sprintf("Channel_Sel_%i", i));
        ch_sel.set("posx", rs.xc(i, 1)).set("posy", rs.xc(i, 2)).set("posz", rs.xc(i, 3));
        ch_sel.set("r", 0);
        ch_sel.set("entitydim", 2);
        channel_indices(i) = ch_sel.entities;

        % flux = tds.feature.create(sprintf("Flux_Ch_%i", i), 'FluxBoundary', 2);
        flux = tds.feature(sprintf("Flux_Ch_%i", i));
        flux.selection.set(channel_indices(i));
        flux.set("species", 1);
        if setting == "complex"
            flux.set('J0', sprintf("-channel_state_%i(t) * P * charge^2 * U * ( (c - (C_extracellular * emU)) / (1 - emU) )", i));
        elseif setting == "block"
            flux.set('J0', sprintf("-ChannelBlock(C_Avg, channel_state_%i(t)) * P * charge^2 * U * ( (c - (C_extracellular * emU)) / (1 - emU) )", i));
        end
    end
end

if setting == "simple"
    channel_open_prob = sum(channel_state_matrix) / length(channel_state_matrix);

    tags = tds.feature.tags;
    for i = 1:length(tags)
        if ~startsWith(tags(i).toCharArray', "Flux_Op_Ch_")
            continue
        end
        tds.feature.remove(tags(i))
    end

    for i = 1:size(rs.xc, 1)
        flux = tds.feature.create(sprintf("Flux_Op_Ch_%i", i), 'FluxBoundary', 2);
        flux.selection.set(channel_indices(i));
        flux.set("species", 1);
        flux.set('J0', sprintf("-%d * P * charge^2 * U * ( (c - (C_extracellular * emU)) / (1 - emU) )", channel_open_prob(i)));
    end
end

if setting == "group"
    % Randomly group 72 channels into 12 groups of 6
    channels_per_group = 6;
    num_groups = 72 / 6;
    channel_indices = reshape(randperm(72), [num_groups, channels_per_group]);

    for group = 1:num_groups
        grouping = channel_indices(group, :);
        group_indices = rs.xc(grouping, :);
        group_sel_name = sprintf('Channel_Group_Sel_%i', group);
        group_flux_name = sprintf('Flux_Op_Group_%i', group);

        % Create selection for the group
        try
            group_sel = model.component(component_name).selection.create(group_sel_name, 'Ball');
        catch
            group_sel = model.component(component_name).selection(group_sel_name);%, 'Ball');
        end
        group_sel.set('entitydim', 2).set("r", 0); % Assuming the channels are defined in 2D

        % Create flux for the group
        try
            flux = tds.feature.create(group_flux_name, 'FluxBoundary', 2);
        catch
            flux = tds.feature(group_flux_name);
        end

        entities = [];
        for i = group_indices'
            group_sel.set('posx', i(1)).set('posy', i(2)).set('posz', i(3));
            entities = [entities group_sel.entities];
        end
        flux.selection.set(entities);


        flux.set('species', 1);
        flux.set('J0', sprintf("-channel_state_%i(t) * P * charge^2 * U * ( (c - (C_extracellular * emU)) / (1 - emU) )", group));
    end
end

end


%% Physics

% channel_open_prob = sum(channel_state_matrix) / length(channel_state_matrix);


% 
% tags = tds.feature.tags;
% for i = 1:length(tags)
%     if ~startsWith(tags(i).toCharArray', "Flux_Op_Ch_")
%         continue
%     end
%     tds.feature.remove(tags(i))
% end

% for i = 1:size(rs.xc, 1)
%     ch_sel = model.component(component_name').selection.create(sprintf("Channel_Sel_%i", i), "Ball");
%     ch_sel.set("posx", rs.xc(i, 1)).set("posy", rs.xc(i, 2)).set("posz", rs.xc(i, 3));
%     ch_sel.set("r", 0);
%     ch_sel.set("entitydim", 2);
%     channel_indices(i) = ch_sel.entities;
% 
%     % flux = tds.feature.create(sprintf("Flux_Op_Ch_%i", i), 'FluxBoundary', 2);
%     % flux.selection.set(channel_indices(i));
%     % flux.set("species", 1);
%     % flux.set('J0', sprintf("-%d * P * charge^2 * U * ( (c - (C_extracellular * emU)) / (1 - emU) )", channel_open_prob(i)));
% end
% 
% tds = model.physics("tds");
% 
% % try
% %     model.component(component_name).variable.create("Average Concentrations");
% % catch
% % end
% % 
% % for i = 1:size(rs.xv, 1)
% %     try
% %         ves_sel = model.component(component_name).selection.create(sprintf("Vesicle_Sel_%i", i), "Ball");
% %     catch
% %         ves_sel = model.component(component_name).selection(sprintf("Vesicle_Sel_%i", i));
% %     end
% %     ves_sel.set("posx", rs.xv(i, 1)).set("posy", rs.xv(i, 2)).set("posz", rs.xv(i, 3));
% %     ves_sel.set("r", 20);
% %     ves_sel.set("entitydim", 2);
% %     % ves_sel.selection.set("condition", "intersects");
% % 
% % 
% %     cpl = model.cpl(sprintf("intop%i", i));
% %     % cpl = model.cpl.create(["intop" i], "Integration");
% %     % cpl.geom("geom", 2);
% %     cpl.selection.set(ves_sel.entities);
% % 
% %     model.component(component_name).variable("Average Concentrations").set(sprintf("C_Avg_%i", i), sprintf("intop%i(c)/intop%i(1)", i, i)).label("Average Concentrations");
% % end
% % try
% %     gev = model.result.numerical.create("gev", "EvalGlobal");
% % catch
% %     gev = model.result.numerical("gev");
% % end
% % for i = 1:size(rs.xv, 1)
% %     try
% %         gev = model.result.numerical.create(sprintf("gev%i", i), "EvalGlobal");
% %     catch
% %         gev = model.result.numerical(sprintf("gev%i", i));
% %     end
% %     gev.set("probetag", "none");
% %     gev.set("expr", sprintf("C_Avg_%i", i));
% %     gev.set("descr", sprintf("Average Calcium Concentration for Vesicle %i", i));
% %     % gev.set("data", "dset2");
% %     % 
% % end
% % gev.setResult;
% 
% 
% 
% if update_interpolation
%     for i = 1:size(rs.xc, 1)
%         func_name = ['channel_state_' num2str(i)];
%         try
%             model.func.create(func_name, 'Interpolation');
%         catch
%         end
% 
%         model.func(func_name).set('interp', 'neighbor');
% 
%         model.func(func_name).set('source', 'table');
% 
%         table_data = cell(length(time_vector), 2);
%         for j = 1:length(time_vector)
%             table_data{j, 1} = num2str(time_vector(j)); % TODO: Move this out of the loop
%             table_data{j, 2} = num2str(channel_state_matrix(j,i)); 
%         end
% 
%         model.func(func_name).set('table', table_data);
% 
%         model.func(func_name).set('argunit', interpolation_unit);
%     end
% end
% 
% 
% % Save the modified model
% try % _MODIFIED
%    model.save(sprintf("C:\\Users\\Tinghan\\Downloads\\ribbon-synapse-main\\%s.mph", file_name));
% catch E
%    disp("Model is locked");
% end
% 
% 
% 


