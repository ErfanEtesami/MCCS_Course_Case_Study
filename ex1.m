%==========================================================================
%   TP :            Case study: Exercse 1
%   Contact:        ezequiel.gonzalezdebada@epfl.ch
%==========================================================================
classdef ex1
    % Class gathering the solutions of exercise 1. 
    methods (Static)
        %
        function varargout = getSystemParameters
            % PARAMETERS = getSystemParameters() returns a 5-elements
            % column vector containing the value of the system parameters 
            % and the linearization point. Specifically it should contain 
            % (in the presented order):
            %   - k : value of curvature of the reference path.
            %   - car_length [m]: car's length.  
            %   - sigma_v : coefficient characterizing the dynamic of the
            %   actuator tracking the speed reference. 
            %   - sigma_phi : coefficient characterizing the dynamic of the
            %   actuator tracking the steering wheel's reference position. 
            %   - spReg : value of the speed the vehicle should drive at. 

            %%- Setting up the parameters' values 
            k = 1e-10; 
            car_length = 4;
            sigma_v = 1;
            sigma_phi = 5;
            spReg = 5;
            %%- Setting up the outputs of the function 
            varargout = {[k; car_length; sigma_v; sigma_phi; spReg]};               
        end
        %
        function varargout = getLinealModelArrays(parameters)
            % [A,B,C,D] = getLinealModelArrays(PARAMETERS) outputs the
            % matrices A,B,C,D characterizing the continuous_time linear
            % version of the system given the set of parameters values as
            % they where defined in function 'getSystemParameters'.
            
            %%- Calculating the matrices A, B, C, D
            k = parameters(1);
            car_length = parameters(2);
            sigma_v = parameters(3);
            sigma_phi = parameters(4);
            spReg = parameters(5);
            A = [0, k*spReg, 0, 1, 0;
                 0, 0, spReg, 0, 0;
                 0, -(k^2)*spReg, 0, 0, (spReg/(16*car_length))*(1+(k^2)*(car_length^2));
                 0, 0, 0, -sigma_v, 0;
                 0, 0, 0, 0, -sigma_phi];
            B = [0, 0;
                 0, 0;
                 0, 0;
                 sigma_v, 0;
                 0, sigma_phi];
            C = eye(5);
            D = zeros(5, 2);
            %%- Setting up the outputs of the function 
            varargout = {[A], [B], [C], [D]};
        end        
        %
        function varargout = getDiscreteLinearModel(A,B,C,D,sampling_time,method)
            % [PHI,GAM] =
            % getDiscreteLinearModel(A, B, C, D, SAMPLING_TIME, METHOD)
            % outputs the PHI and GAMMA matrices characterizing the
            % discrete linear model of the system given the matrices
            % A,B,C,D of the continuous time linear description of the
            % system and the desired SAMPLING_TIME.
            %
            % Additionally, the input METHOD will contain a string
            % indicating the way the matrices PHI and GAMMA are wanted to
            % be calculated. Such an input can take values 
            % - Euler : Euler approximation as discretization method. 
            % - Psi : Use the progrmmatically algorithm presented in the
            % note course which makes use of the intermediary matrix Psi.
            % - c2d : use the matlab command c2d. 
            
            %%- Implementing discretization methods based on the inputs
            if strcmp(method,'Euler')
                % Use here Euler approximation to approximate the
                % discrete-time linear model of the system.
                Phi = eye(size(A, 1)) + sampling_time * A;
                Gamma = sampling_time * B;
                %%- Setting up the outputs of the function
                varargout = {[Phi], [Gamma]};
            elseif strcmp(method,'Psi')
                % Dimension of Psi and Psi initialization
                psi_dimention = size(A, 1); 
                Psi = zeros(psi_dimention);
                number_of_iterations = 10;  % Arbitrarily chosen
                for i = 0:1:number_of_iterations
                    % Updating Psi iteratively 
                    Psi = Psi + (A^i)*(sampling_time^i)/factorial(i+1);
                end
                % Calculateing matrices Phi and Gamma
                Phi = eye(size(A, 1)) + sampling_time * A * Psi;
                Gamma = sampling_time * Psi * B;
                %%- Setting up the outputs of the function
                varargout = {[Phi], [Gamma]};                         
            elseif strcmp(method,'c2d')
                % Continuous representation of the system with 'ss'
                Mc = ss(A, B, C, D);
                % Calculateing the discrete-time linear model of the system
                % using 'c2d'
                Md = c2d(Mc, sampling_time); 
                % Extracting the Phi and Gamma matrices from Md 
                Phi = Md.A;
                Gamma = Md.B;
                %%- Setting up the outputs of the function
                varargout = {[Phi], [Gamma]};                
            end
        end                
        %
        function varargout = getWorkingTrajectory(sampling_time, simulation_time, parameters)
            % [NOMINAL_TRAJECTORY_X, NOMINAL_TRAJECTORY_U] =
            % getWorkingTrajectory(SAMPLING_TIME, SIMULTAION_TIME,
            % PARAMETERS)  
            % outputs the NOMINAL_TRAJECTORY_X and NOMINAL_TRAJECTORY U
            % given the SAMPLING_TIME between data points, the
            % SIMULATION_TIME up to which the trajectory has to be created,
            % and the vector PARAMETERS with the value sof tha system's
            % parameters.
            %
            % Outputs NOMINAL_TRAJECTORY_X, and NOMINAL_TRAJECTORY_U must
            % be arrays whose first collumn correspond to the time span of
            % the data point, and successive columns store the information
            % of the states and inputs, correspondingly.
            %
            % The defined output trajectories are meant to be used in
            % Simulink with "From Workspace" importing module. If any
            % additional doubt regarding how the data should be structured,
            % read the information provided by the mentioned simulink block.
            %
            % To Do:
            % - create time vector. 
            % - create the nominal states trajectory output
            % - create the control inputs nominal trajectory output
            
            %%- Creating time vector
            time_vector = [0:sampling_time:simulation_time]';
            %%- Creating nominal state trajectory
            x1_bar_dot = parameters(5);
            x1_bar = x1_bar_dot * time_vector;
            x2_bar = zeros(length(time_vector), 1);
            x3_bar = zeros(length(time_vector), 1);
            x4_bar = parameters(5) * ones(length(time_vector), 1);
            x5_bar = 16*atan(parameters(1)*parameters(2)) * ones(length(time_vector), 1);
            nominal_trajectory_x = [time_vector, x1_bar, x2_bar, x3_bar, x4_bar, x5_bar];
            %%- Creating nominal control input trajectory 
            u1_bar = x4_bar;
            u2_bar = x5_bar;
            nominal_trajectory_u = [time_vector, u1_bar, u2_bar];
            %%- Setting up the outputs of the function 
            varargout = {[nominal_trajectory_x], [nominal_trajectory_u]};
        end
        %
        function varargout = getInitialState(nominal_trajectory_x)
            %[X0, X0TILDE] = getInitialState(NOMINAL_TRAJECTORY_X)
            % returns the initial state X0 of the non linear system as the
            % initial state X0TILDE of the linearized model of the system,
            % given the information on the exercise handout and the
            % NOMINAL_TRAJECTORY_X.
            %
            % The outputs should be column vectors. 
            %
            % Remember that by definition \tilde{x} = x - \overline{x}.
            
            %%- Defining the value of x0 for experiment 1
            x0_experiment_1 = [nominal_trajectory_x(1, 2:end)]';
            % x0_experiment_1 = [0; 0; 0; nominal_trajectory_x(1, 5); nominal_trajectory_x(1, 6)];
            %%- Defining the value of x0Tilde for experiment 1
            x0Tilde_experiment_1 = x0_experiment_1 - [nominal_trajectory_x(1, 2:end)]'; 
            %%- Including the different values for different experiments 
            %%- as a cell
            x0 = {x0_experiment_1};
            x0Tilde = {x0Tilde_experiment_1};
            %%- Setting up the outputs of the function 
            varargout = {[x0],[x0Tilde]};
        end
        %
        function varargout = getOpenLoopControlSignal(sampling_time, simulation_time)
            %[INPUT_CONTROL_ACTIONS_OPEN_LOOP] = getOpenLoopControlSignal(SAMPLING_TIME, SIMULATION_TIME)
            % outputs a sequence of control signals to be applied in open
            % loop to the system, which will be then simulated. In order to
            % do that, the desired SAMPLING_TIME between data points as
            % well as the SIMULTION_TIME is provided. 
            %
            % As int he case of GETWORKINGTRAJECTORY function, the outputs
            % are meant to be used in Simulink with "From Workspace"
            % importing module. If any additional doubt regarding how the
            % data should be structured, read the information provuded by
            % the mentioned simulink block. 
            %
            % To Do:
            % - Declare an appropriate time span vector. 
            % - Create the input_control_actions_open_loop array with the
            % sequence of control inputs to be applied in open loop. 
            %
            % Notice: alternatively, this function can output a cell with
            % several arrays showing different control sequences to be
            % applied. This would make the provided script to run the
            % simulink model as many times as different control sequences
            % are gathered within the cell. Meaning that several
            % experiments can be set at once. 
            
            %%- Creating a time vector
            time_vector = [0:sampling_time:simulation_time]';
            %%- Setting the control sequence to be applied in open loop for 
            %%- the 1st experiment
            v_ref = 5;
            u1 = v_ref * ones(length(time_vector), 1);
            phi_ref = 0;
            u2 = phi_ref * ones(length(time_vector), 1);
            uOpenLoop_experiment = [time_vector, u1, u2];
            %%- the 2nd experiment
%             v_ref = 5;
%             u1 = v_ref * time_vector;
%             phi_ref = 0;
%             u2 = phi_ref * ones(length(time_vector), 1);
%             uOpenLoop_experiment = [time_vector, u1, u2];
            %%- the 3rd experiment
%             v_ref = 5;
%             u1 = v_ref * time_vector;
%             phi_ref = pi / 36;
%             u2 = phi_ref * ones(length(time_vector), 1);
%             uOpenLoop_experiment = [time_vector, u1, u2];
            %%- the 4th experiment
%             v_ref = 5;
%             u1 = v_ref * time_vector;
%             phi_ref = pi / 6;
%             u2 = phi_ref * ones(length(time_vector), 1);
%             uOpenLoop_experiment = [time_vector, u1, u2];
            %%- the 5th experiment
%             v_ref = 5;
%             u1 = v_ref * time_vector;
%             phi_ref = 2 * pi;
%             u2 = [linspace(0, phi_ref, length(time_vector))]';
%             uOpenLoop_experiment = [time_vector, u1, u2];
            %%- the 6th experiment
%             v_ref = 5;
%             u1 = [linspace(0, v_ref, length(time_vector))]';
%             phi_ref = 2 * pi;
%             u2 = [linspace(0, phi_ref, length(time_vector))]';
%             uOpenLoop_experiment = [time_vector, u1, u2];
            %%- Including different values for the different experiments in 
            %%- a cell
            input_control_actions_open_loop = {uOpenLoop_experiment};
            %%- Setting up the output of the function
            varargout = {[input_control_actions_open_loop]};
        end
        % 
    end
    %
end
%EOF