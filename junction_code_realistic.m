close all
clearvars

tic

%This code randomly chooses connection point amongst all points within the
%given threshold, weighted by inverse distance (more likely to connect at
%the closest points). Diffusion of rods from Broersma theory is used for
%translational and rotational diffusion.

%Last updated: 7/16/25

%% Parameters
L = 115; % side length of the cube in um
segment_length = 12; % length of each line segment in um
num_segments = 200; % number of line segments
Tend = 420; %length of time in seconds to run simulation for
dt = 0.01; % time step in s
num_steps = ceil(Tend/dt); % number of steps in the simulation
%rotational_diffusion = 0.2; % scale factor of the random orientation (rotational) movement
%translational_diffusion = 1; %scale factor of the random translational movement
proximity_threshold = 0.1; % threshold distance for segments to stick together in um
discrete_size = 20; % number of discretization points per segment
plot_interval = 1000;

%Parameters for diffusion coefficients
kB = 1.380649*10^(-23); %Boltzmann's constant (m^2 * kg) / (s^2 * K)
T = 310.15; %human body temperature in K (room temp, which was experiments, is about 296)
eta = 1.8*6.193*10^(-4); %dynamic viscosity of blood plasma is 1.8-times the viscosity of water in kg / (m * s)
drod = 1*10^(-7); %diameter of rod in m (=100 nm)
Lrod = segment_length*10^(-6); %length of rod in m

%calculate translational diffusion coefficient - Broersma corrected
Dtrans = kB*T*(log(Lrod/drod)+(0.866-0.15/(log(2*Lrod/drod))-8.1/(log(2*Lrod/drod))^2+18/(log(2*Lrod/drod))^3-9/(log(2*Lrod/drod))^4-0.114-0.15/(log(2*Lrod/drod))-13.5/(log(2*Lrod/drod))^2+37/(log(2*Lrod/drod))^3-22/(log(2*Lrod/drod))^4)/2)/(3*pi*eta*Lrod); %measured in m^2/s
%Dtrans = kB*T*log(Lrod/drod)/(6*pi*eta*Lrod); %measured in m^2/s
Dtrans = Dtrans*10^(12); %measured in um^2/s

%calculate rotational diffusion coefficient - Broersma corrected
Drot = 3*kB*T*(log(Lrod/drod)-0.446-0.2/log(2*Lrod/drod)-16/(log(2*Lrod/drod))^2+63/(log(2*Lrod/drod))^3-62/(log(2*Lrod/drod))^4)/(pi*eta*Lrod^3); %measured in rad^2/s
%Drot = 3*kB*T/(pi*eta*Lrod^3); %measured in rad^2/s

%rng(1) % Gives the same random numbers every time, comment out if needed

%% Initialize the positions of the endpoints of the line segments
% Generate random centers for segments with margin for segment length
margin = segment_length / 2;
centers = rand(num_segments, 3) * (L - 2*margin) + margin;

% Generate random directions
directions = randn(num_segments, 3); % random initial directions
directions = directions ./ vecnorm(directions, 2, 2); % normalize directions

% Calculate endpoints from centers
p1 = centers - (segment_length/2) * directions;
p2 = centers + (segment_length/2) * directions;

% Initialize cell array for discretized segments
P = cell(1, num_segments);
% Discretize each segment
for i = 1:num_segments
    P{i} = [linspace(p1(i,1), p2(i,1), discrete_size);
        linspace(p1(i,2), p2(i,2), discrete_size);
        linspace(p1(i,3), p2(i,3), discrete_size)];
end

% Initialize connection tracking - stores [seg1, seg2, point1, point2] for each connection
connections = zeros(num_segments^2, 4);
connections_len = 0;
junction_type = zeros(num_segments^2);
junction_type_len = 0;
theta = zeros(num_segments^2);
theta_len = 0;

% Initialize disjoint-set data structure
[Pred Rank] = InitializeRootedTree(num_segments)
perpetual_components = cell(num_segments, 1)
for i = 1:num_segments
    perpetual_components{i} = [i]
end

%% Plotting setup
figure(1);
view(3);
hold on;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([0 L]); ylim([0 L]); zlim([0 L]);
title('Fiber Diffusion and Junction Formation in 3D');
colors = lines(num_segments);

% Initialize plot handles
% segment_plots = zeros(num_segments, 1);
% for i = 1:num_segments
%     segment_plots(i) = plot3(P{i}(1,:), P{i}(2,:), P{i}(3,:), '-', ...
%         'Color', colors(i,:), 'LineWidth', 2);
% end
% drawnow;



%% Main simulation loop
for step = 1:num_steps

    % Move each component as a rigid body
    for comp_idx = 1:length(perpetual_components)
        component = perpetual_components{comp_idx};

        if length(component) == 0
            continue;
        end

        % Generate random movement and rotation for this component
        noise_trans = randn(1,3)*sqrt(2*Dtrans*dt);
        noise_rot = randn(1,3)*sqrt(2*Drot*dt);
        

        if length(component) == 1
            % Single segment - move normally
            i = component(1);
            center = mean(P{i}, 2);
            current_direction = P{i}(:,end) - P{i}(:,1);
            current_direction = current_direction / norm(current_direction);

            % Apply random movement to center
            new_center = center + noise_trans';

            % Apply random orientation change
            axis = noise_rot';
            angle = norm(axis);
            if angle > 0
                axis = axis / angle;
                new_direction = current_direction + angle * cross(axis, current_direction);
                new_direction = new_direction / norm(new_direction);
            else
                new_direction = current_direction;
            end

            % Update segment endpoints
            p1_new = new_center - (segment_length/2) * new_direction;
            p2_new = new_center + (segment_length/2) * new_direction;

            % Apply boundary conditions
            p1_new = max(0, min(L, p1_new));
            p2_new = max(0, min(L, p2_new));

            % Ensure segment length is maintained
            actual_direction = p2_new - p1_new;
            actual_length = norm(actual_direction);
            if actual_length > 0
                actual_direction = actual_direction / actual_length;
                p2_new = p1_new + segment_length * actual_direction;
            end

            % Re-discretize the segment
            P{i} = [linspace(p1_new(1), p2_new(1), discrete_size);
                linspace(p1_new(2), p2_new(2), discrete_size);
                linspace(p1_new(3), p2_new(3), discrete_size)];
        else
            % Connected component - move as rigid body while preserving connections

            % Store original connection points before any movement
            original_connection_points = [];
            for conn_idx = 1:connections_len
                seg1 = connections(conn_idx, 1);
                seg2 = connections(conn_idx, 2);
                point1_idx = connections(conn_idx, 3);
                point2_idx = connections(conn_idx, 4);

                % Only store connections within this component
                if ismember(seg1, component) && ismember(seg2, component)
                    original_connection_points(end+1,:) = [conn_idx, seg1, seg2, point1_idx, point2_idx];
                end
            end

            % Calculate center of mass of the component
            all_points = [];
            for i = component
                all_points = [all_points, P{i}];
            end
            com = mean(all_points, 2);

            % Apply translation
            % new_com = com + (1/length(perpetual_components{comp_idx}))*noise_trans'; %have bigger components move slower
            % for i = component
            %     P{i} = P{i} + (1/length(perpetual_components{comp_idx}))*noise_trans'; %have bigger components move slower
            % end
            new_com = com + (1/(length(perpetual_components{comp_idx}))^3)*noise_trans'; %have bigger components move MUCH slower
            for i = component
                P{i} = P{i} + (1/(length(perpetual_components{comp_idx}))^3)*noise_trans'; %have bigger components move MUCH slower
            end

            % Apply rotation around center of mass
            %axis = (1/length(perpetual_components{comp_idx}))*noise_rot'; %have bigger components move slower
            axis = (1/(length(perpetual_components{comp_idx}))^3)*noise_rot'; %have bigger components move MUCH slower
            angle = norm(axis);
            if angle > 0
                axis = axis / angle;

                % Create rotation matrix for small angles
                K = [0, -axis(3), axis(2);
                    axis(3), 0, -axis(1);
                    -axis(2), axis(1), 0];
                R = eye(3) + sin(angle) * K + (1-cos(angle)) * K^2;

                % Apply rotation to each segment around component
                for i = component
                    P{i} = R * (P{i} - new_com) + new_com;
                end
            end

            % Enforce connection constraints after movement
            for conn_row = 1:size(original_connection_points, 1)
                conn_idx = original_connection_points(conn_row, 1);
                seg1 = original_connection_points(conn_row, 2);
                seg2 = original_connection_points(conn_row, 3);
                point1_idx = original_connection_points(conn_row, 4);
                point2_idx = original_connection_points(conn_row, 5);

                % Get current positions of connection points
                point1 = P{seg1}(:, point1_idx);
                point2 = P{seg2}(:, point2_idx);

                % Calculate midpoint
                midpoint = (point1 + point2) / 2;

                % Move both connection points to the midpoint
                displacement1 = midpoint - point1;
                displacement2 = midpoint - point2;

                % Apply constraint forces to maintain segment lengths
                % For segment 1
                p1_seg1 = P{seg1}(:, 1);
                p2_seg1 = P{seg1}(:, end);
                center1 = (p1_seg1 + p2_seg1) / 2;

                % Move segment 1 to satisfy connection constraint
                P{seg1} = P{seg1} + displacement1;

                % Ensure segment 1 maintains its length
                new_p1_seg1 = P{seg1}(:, 1);
                new_p2_seg1 = P{seg1}(:, end);
                current_dir1 = new_p2_seg1 - new_p1_seg1;
                current_length1 = norm(current_dir1);
                if current_length1 > 0
                    current_dir1 = current_dir1 / current_length1;
                    new_center1 = (new_p1_seg1 + new_p2_seg1) / 2;
                    new_p1_seg1 = new_center1 - (segment_length/2) * current_dir1;
                    new_p2_seg1 = new_center1 + (segment_length/2) * current_dir1;

                    % Re-discretize segment 1
                    P{seg1} = [linspace(new_p1_seg1(1), new_p2_seg1(1), discrete_size);
                        linspace(new_p1_seg1(2), new_p2_seg1(2), discrete_size);
                        linspace(new_p1_seg1(3), new_p2_seg1(3), discrete_size)];
                end

                % For segment 2
                p1_seg2 = P{seg2}(:, 1);
                p2_seg2 = P{seg2}(:, end);
                center2 = (p1_seg2 + p2_seg2) / 2;

                % Move segment 2 to satisfy connection constraint
                P{seg2} = P{seg2} + displacement2;

                % Ensure segment 2 maintains its length
                new_p1_seg2 = P{seg2}(:, 1);
                new_p2_seg2 = P{seg2}(:, end);
                current_dir2 = new_p2_seg2 - new_p1_seg2;
                current_length2 = norm(current_dir2);
                if current_length2 > 0
                    current_dir2 = current_dir2 / current_length2;
                    new_center2 = (new_p1_seg2 + new_p2_seg2) / 2;
                    new_p1_seg2 = new_center2 - (segment_length/2) * current_dir2;
                    new_p2_seg2 = new_center2 + (segment_length/2) * current_dir2;

                    % Re-discretize segment 2
                    P{seg2} = [linspace(new_p1_seg2(1), new_p2_seg2(1), discrete_size);
                        linspace(new_p1_seg2(2), new_p2_seg2(2), discrete_size);
                        linspace(new_p1_seg2(3), new_p2_seg2(3), discrete_size)];
                end
            end

            % Apply boundary conditions to entire component
            for i = component
                % Get current endpoints
                p1_current = P{i}(:,1);
                p2_current = P{i}(:,end);

                % Apply boundary reflection to endpoints
                p1_new = p1_current;
                p2_new = p2_current;

                % Reflect if outside boundaries
                for dim = 1:3
                    if p1_current(dim) < 0
                        p1_new(dim) = -p1_current(dim);
                    elseif p1_current(dim) > L
                        p1_new(dim) = 2*L - p1_current(dim);
                    end

                    if p2_current(dim) < 0
                        p2_new(dim) = -p2_current(dim);
                    elseif p2_current(dim) > L
                        p2_new(dim) = 2*L - p2_current(dim);
                    end
                end

                % Maintain segment length
                direction = p2_new - p1_new;
                current_length = norm(direction);
                if current_length > 0
                    direction = direction / current_length;
                    p2_new = p1_new + segment_length * direction;
                end

                % Re-discretize to maintain straight line
                P{i} = [linspace(p1_new(1), p2_new(1), discrete_size);
                    linspace(p1_new(2), p2_new(2), discrete_size);
                    linspace(p1_new(3), p2_new(3), discrete_size)];
            end
        end
    end

    % Check for proximity and create NEW connections (only for unconnected segments)
    for i = 1:num_segments
        for j = i+1:num_segments
            root_i = FindRoot(Pred, i)
            root_j = FindRoot(pred, j)

            if root_i == root_j
                continue;
            end

            % Calculate all distances between segments and find candidates within threshold
            candidates_len = 0;
            candidates = zeros(num_segments^2, 3);
            for pI = 1:discrete_size
                for pj = 1:discrete_size
                    dist = norm(P{i}(:,pI) - P{j}(:,pj));
                    if dist < proximity_threshold
                        candidates_len++
                        candidates(candidates_len, :) = [pI, pj, dist];
                    end
                end
            end

            % If we have candidates, select one using distance-weighted probability
            if candidates_len > 0
                candidates = candidates(1:candidates_len, :)
                % Calculate weights (inverse distance, so closer pairs are more likely)
                weights = 1 ./ (candidates(:,3) + 1e-10); % Add small epsilon to avoid division by zero
                weights = weights / sum(weights); % Normalize to probabilities

                % Select a candidate based on weighted probability
                cumulative_weights = cumsum(weights);
                rand_val = rand();
                selected_idx = find(cumulative_weights >= rand_val, 1);

                closest_point_i = candidates(selected_idx, 1);
                closest_point_j = candidates(selected_idx, 2);
                min_dist = candidates(selected_idx, 3);

                % Calculate midpoint between the selected points
                point_i = P{i}(:, closest_point_i);
                point_j = P{j}(:, closest_point_j);
                midpoint = (point_i + point_j) / 2;

                % Move both segments so their connection points meet at the midpoint
                displacement_i = midpoint - point_i;
                displacement_j = midpoint - point_j;

                P{i} = P{i} + displacement_i;
                P{j} = P{j} + displacement_j;

                % Union the two sets
                [Pred Rank] = UnionbyRank(Pred, Rank, root_i, root_j)
                if root_i == FindRoot(Pred, root_i)
                    perpetual_components{root_i} = [perpetual_components{root_i} perpetual_components{root_j}]
                    perpetual_components{root_j} = []
                else
                    perpetual_components{root_j} = [perpetual_components{root_j} perpetual_components{root_i}]
                    perpetual_components{root_i} = []
                end

                % Add connection to matrix
                connections_len++
                connections(connections_len, :) = [i, j, closest_point_i, closest_point_j];

                % Calculate angle between newly connected fibers
                new_angle = acos(dot((P{i}(:,discrete_size)-P{i}(:,1)),P{j}(:,discrete_size)-P{j}(:,1))/(norm(P{i}(:,discrete_size)-P{i}(:,1))*norm(P{j}(:,discrete_size)-P{j}(:,1))));
                theta_len++
                theta(theta_len) = new_angle; %measured in radians

                % Determine connection type for reporting
                is_endpoint_i = (closest_point_i == 1) || (closest_point_i == discrete_size);
                is_endpoint_j = (closest_point_j == 1) || (closest_point_j == discrete_size);

                connection_type = '';
                if is_endpoint_i && is_endpoint_j
                    connection_type = ' (endpoint-endpoint)';
                    new_junc = 1; %1 means end-to-end, or L-junction
                elseif is_endpoint_i || is_endpoint_j
                    connection_type = ' (endpoint-interior)';
                    new_junc = 2; %2 means end-to-middle, or T-junction
                else
                    connection_type = ' (interior-interior)';
                    new_junc = 3; %3 means middle-to-middle, or X-junction
                end

                fprintf('Step %d: Segments %d and %d connected at points %d and %d%s (distance: %.3f -> 0.000) (angle: %.3f)\n', ...
                    step, i, j, closest_point_i, closest_point_j, connection_type, min_dist, new_angle);
                junction_type_len++
                junction_type(junction_type_len)=new_junc;
            end
        end
    end

    % Update plots
    % if mod(step,plot_interval) == 0 %only plot every "plot_interval" time steps
    %     step
    %     for i = 1:num_segments
    %         set(segment_plots(i), 'XData', P{i}(1,:), ...
    %             'YData', P{i}(2,:), ...
    %             'ZData', P{i}(3,:));
    %     end

    %     % Visualize connections - draw both connection points
    %     for conn_idx = 1:connections_len
    %         seg1 = connections(conn_idx, 1);
    %         seg2 = connections(conn_idx, 2);
    %         point1 = connections(conn_idx, 3);
    %         point2 = connections(conn_idx, 4);

    %         % Verify that connection points are still coincident
    %         pos1 = P{seg1}(:, point1);
    %         pos2 = P{seg2}(:, point2);
    %         connection_error = norm(pos1 - pos2);

    %         if connection_error > 1e-3 % Small tolerance for numerical errors
    %             fprintf('Warning: Connection %d has error %.6f\n', conn_idx, connection_error);
    %         end

    %         % Draw connection points
    %         plot3(pos1(1), pos1(2), pos1(3), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
    %         plot3(pos2(1), pos2(2), pos2(3), 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'green');
    %     end

    %     drawnow;

    %     % Capture frame
    %     % if step == 1
    %     %     frames = getframe(gcf);
    %     % else
    %     %     frames(step) = getframe(gcf);
    %     % end

    %     % Clear connection markers for next frame
    %     children = get(gca, 'Children');
    %     markers_to_delete = [];
    %     for c = 1:length(children)
    %         if strcmp(get(children(c), 'Marker'), 'o')
    %             markers_to_delete = [markers_to_delete, c];
    %         end
    %     end
    %     delete(children(markers_to_delete));
    % end
end

toc

% Display final statistics
total_connections = connections_len;

fprintf('\nSimulation completed!\n');
fprintf('Total connections formed: %d\n', total_connections);
%fprintf('Frames captured: %d\n', length(frames));

% Trim off the excess of the junction_type array
junction_type = junction_type(1:junction_type_len)

fprintf('number of X-junctions: %d\n',length(find(junction_type==3)));
fprintf('number of T-junctions: %d\n',length(find(junction_type==2)));
fprintf('number of L-junctions: %d\n',length(find(junction_type==1)));

% Trim off the excess of the theta array
theta = theta(1:theta_len)
%convert angle from radians to degrees
thetadeg=theta*180/pi;

%only use the acute angles
thetadegacute=min(thetadeg, 180 - thetadeg);

components = {}
for i = 1:length(perpetual_components)
    if length(perpetual_components{i}) > 0
        components{end+1} = perpetual_components{i}
    end
end

%save the data from the run
save -ascii realistic_Tend420_200fibers_slow_diffusion_junction_type.dat junction_type
save -ascii realistic_Tend420_200fibers_slow_diffusion_junction_angle.dat thetadegacute
save realistic_Tend420_200fibers_slow_diffusion_components.mat components
save('realistic_Tend420_200fibers_slow_diffusion_parameters.mat', 'L', 'segment_length','num_segments','Tend','dt','proximity_threshold','discrete_size','Dtrans','Drot');


%visualize histograms of angles
figure(2)
histogram(thetadegacute,'BinWidth',5)
title('Histogram of angles')
savefig('angles_realistic_Tend420_200fibers_slow_diffusion.fig')

figure(3)
histogram(junction_type)
title('Histogram of junction types')
savefig('junction_type_realistic_Tend420_200fibers_slow_diffusion.fig')




function [Pred Rank] = UnionbyRank(Pred, Rank, Tree1, Tree2)
% This function accepts a Predecessor and Rank function describing a rooted
% tree forest and two trees to be merged. It then merges the trees based on
% their ranks and updates the Pred and Rank lists

if Rank(Tree1) > Rank(Tree2)
    Pred(Tree2) = Tree1; % Merge the smaller tree into the larger tree
elseif Rank(Tree1) < Rank(Tree2)
    Pred(Tree1) = Tree2; % Merge the smaller tree into the larger tree
else
    Pred(Tree1) = Tree2; % Merge the trees arbitrarily
    Rank(Tree2) = Rank(Tree2) + 1; % Increase the rank of the new tree by 1
end

function [Pred Root] = FindRoot(Pred, Node)
% This function accepts A Predecessor list for a rooted forest and a Node.
% The function then finds the root of the tree containing "Node" and, in
% the process, uses path halving to shorten subsequent searches.

Root = Node; % Set the root equal to the start Node
while Pred(Root) ~= Root % as long as the predecessor of the root is not the root itself
    OldNode = Root; % Store the Node we're currently at
    Root = Pred(Root); % Move one step up the tree
    Pred(OldNode) = Pred(Root); % Set the predecessor of the stored node to the node two levels above it
    Root = Pred(Root); % Move another step up the tree
end
Root = Pred(Root); % Once the "while" loop exits we could still be one level below the root so we'll move up one more step.
% NOTE: we may make up to two more steps than we need to but this won't
% be a problem since Pred(Root) = Root so we'll just cycle for a few steps.
% Also the additional steps caused by this are well worth the savings from
% the path halving.

function [Pred Rank] = InitializeRootedTree(n)
% This function creates a predecessor list and a rank list for "n" nodes in
% a rooted tree forest. At this point all nodes are in their own tree so
% pred(i)=i and Rank(i)=0

for i = 1:n
    Pred(i) = i;
    Rank(i) = 0;
end

