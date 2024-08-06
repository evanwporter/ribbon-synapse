rng(21, "twister")

directories = dir('Data/*');
directories = directories(5:end);

total_release = [];
folders = {};

for k = 1:length(directories)
    folder = directories(k).name;
    folders{k} = folder;

    release = [];
    voltage_steps_actual = [];

    for Vt = voltage_steps 
    % for Vt = voltage_steps
        Vt_STR = num2str(abs(Vt));

        try
            tbl_data = readmatrix(sprintf("Data\\%s\\Vt_%s.csv", folder, Vt_STR(3:end)));
            voltage_steps_actual(end + 1) = Vt;
        catch
            continue
        end
        % Vt = -0.07;
        % tbl_data = readmatrix("Data\\Vt_004.csv");

        conc_matrix = tbl_data(:, 3:end);
        num_vesicles = size(conc_matrix, 2);

        q = ones(num_vesicles, 1);
        % q(6) = 0;
        c = 0;
        w = 0;

        dt = 1e-4;

        q_array = ones(num_vesicles, length(conc_matrix)+1);
        amount_released = zeros(length(conc_matrix), 1);
        m = 0;

        for i = 1:length(conc_matrix)
            [ dq, dc, dw, released] = NTdynamics( q, c, w, num_vesicles, opts, dt, conc_matrix(i, :));
            q = q + dq * dt;
            c = c + dc * dt;
            w = w + dw * dt;

            amount_released(i) = sum(released);

            q_array(:, i+1) = q;

            % if any(dq)
            %     disp("NON ZERO")
            % end
        end

        disp(Vt);
        disp(["RELEASED PER SECOND: " num2str(sum(amount_released))]);

        release = [release, sum(amount_released)];

        % release_T{k} = amount_released;

    end

    % total_release = [total_release; release];

    data{k} = [voltage_steps_actual; release];

    plot(data{k});
    hold on;
end