generate_coordinates(2)
generate_coordinates(7)

function generate_coordinates(tightness)
    % Number of objects
    num_objects = 72;
    half_num_objects = num_objects / 4; % Number of objects per quadrant

    % Define the oval dimensions based on tightness
    max_radius_x = 300; % Max x-axis radius
    max_radius_y = 150;  % Max y-axis radius

    % Calculate the effective radii based on tightness
    radius_x = max_radius_x * tightness;
    radius_y = max_radius_y * tightness;

    % Generate symmetrical coordinates within the oval shape
    coords = zeros(num_objects, 2);
    for i = 1:half_num_objects
        while true
            % Generate random coordinates within the bounding box in the first quadrant
            x = rand * radius_x;
            y = rand * radius_y;
            
            % Check if the point is inside the oval
            if (x / radius_x)^2 + (y / radius_y)^2 <= 1
                coords(i, :) = [x, y];
                coords(half_num_objects + i, :) = [-x, y]; % Reflect across the y-axis
                coords(2 * half_num_objects + i, :) = [x, -y]; % Reflect across the x-axis
                coords(3 * half_num_objects + i, :) = [-x, -y]; % Reflect across both axes
                break;
            end
        end
    end

    % Adjust coordinates to fit with blue spheres
    coords(:, 1) = coords(:, 1) + 300; % Centering the distribution
    coords(:, 2) = coords(:, 2) + 50; % Centering the distribution

    % Blue sphere coordinates
    blue_spheres = [
         0   100    60;
         0     0    60;
       100   100    60;
       100     0    60;
       200   100    60;
       200     0    60;
       300   100    60;
       300     0    60;
       400   100    60;
       400     0    60;
       500   100    60;
       500     0    60;
       600   100    60;
       600     0    60;
       700     0    60;
       700   100    60
    ];

    % Plot the results for visualization
    figure;
    hold on;

    % Plot the oval for reference
    t = linspace(0, 2*pi, 100);
    x_oval = radius_x * cos(t) + 350; % Center the oval along the x-axis
    y_oval = radius_y * sin(t) + 50;
    plot(x_oval, y_oval, 'k--'); % Oval boundary

    % Plot the red cylinders and blue spheres
    scatter(coords(:, 1), coords(:, 2), 'filled', 'r'); % Red cylinders
    scatter(blue_spheres(:, 1), blue_spheres(:, 2), 'filled', 'b'); % Blue spheres

    axis equal;
    title(['Object Distribution with Tightness Parameter = ', num2str(tightness)]);
    xlabel('X');
    ylabel('Y');
    legend({'Oval Boundary', 'Red Cylinders', 'Blue Spheres'});
    grid on;
    hold off;
end
