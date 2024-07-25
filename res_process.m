folderPath = 'Results';
matFiles = dir(fullfile(folderPath, '*.mat'));
NTreleased = [];

% Loop through each .mat file
for k = 1:length(matFiles)
    matFileName = fullfile(folderPath, matFiles(k).name);
    loadedData = load(matFileName);

    if isfield(loadedData, 'y_out')
        data = loadedData.data;
        result = calc_q_released(data);
        NTreleased = [NTreleased; result];
        clear loadedData data;
    else
        disp(['Variable "data" not found in ', matFiles(k).name]);
    end
end

figure;
plot(voltage_steps, release_rates, '-o');
xlabel('Voltage (V)');
ylabel('Neurotransmitter Release Rate');
title('Neurotransmitter Release Rate vs Voltage');
grid on;

function total_release = calc_q_released(t, z_array, opts)
    q_start = opts.size_info.NT_free.start;
    q_end = opts.size_info.NT_free.end;

    q_values = z_array(:, q_start:q_end);

    % Calculate total release by summing the decreases in q
    total_release = 0;
    for i = 2:length(t)
        try
            released_this_step = q_values(i - 1, :) - q_values(i, :);
        catch e
            throw(e)
        end
        total_release = total_release + sum(released_this_step(released_this_step > 0));
    end

    disp("Total neurotransmitter released: " + total_release);
end