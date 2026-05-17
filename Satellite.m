classdef Satellite < handle
    % Class for forming objects of type satellite, managing their orbital
    % elements/state, determining missile observability, and forming
    % measurements and simulated measurements.
    
    properties
        orbElems % [a; e; i; O; w; f] initially, after that the state "x" is used for everything
        x % [px; py; pz; vx; vy; vz] in ECI

        constellationName % string denoting which constellation satellite belongs to
        satID % int denoting the satellite ID

        z % IR sensor measurement
        y % predicted IR sensor measurement
        canObserve = false; % boolean denoting whether the satellite can observe the ballistic missile at any given time
        H_func
        H

        hasLOS = false;
        meetsIntensity = false;
        IR_sensor_gain = 1; % make IR sensors more or less sensitive to detection
        I_threshold = 1e-3; % TODO edit
        

        mu = 3.986e14; % Earth's gravitational constant, m^3/s^2
    end
    
    methods
        function obj = Satellite(satID, constellationName, a, e, i, O, w, f)

            obj.satID = satID;
            obj.constellationName = constellationName;

            obj.orbElems = [a;e;i;O;w;f];
            [r, v] = obj.elements_to_rv(a, e, i, O, w, f, obj.mu);
            obj.x = [r;v];

            obj.formHfunction();

            
        end

        function [xdot]= satelliteDynamics(obj, x)

            r = x(1:3);
            v = x(4:6);
            xdot = zeros(6,1);
            xdot(1:3) = v;
            xdot(4:6) = (-obj.mu/(norm(r)^3))*r;

        end

        function [obj]= propagateSatellite(obj, dt)
            % Updates object state obj.x with satellite dynamics function
            % and RK4 (faster in a large scale simulation than ode45())
                
            k1 = obj.satelliteDynamics(obj.x);
            k2 = obj.satelliteDynamics(obj.x + (dt/2)*k1);
            k3 = obj.satelliteDynamics(obj.x + (dt/2)*k2);
            k4 = obj.satelliteDynamics(obj.x + dt*k3);

            obj.x = obj.x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        end

        function [obj]= checkLOS(obj, r_missile_ECI)
            % Function to check if missile has line of sight to satellite

            r_satellite_ECI = obj.x(1:3);

            lambda = 0:0.01:1;

            line_to_missile = r_satellite_ECI + (r_missile_ECI-r_satellite_ECI)*lambda;

            norm_line = vecnorm(line_to_missile, 2);
            min_altitude_visible = 6371000 + 20000; % Radius of Earth plus where atmosphere is thick

            if any(norm_line <= min_altitude_visible)
                obj.hasLOS = false;
            else
                obj.hasLOS = true;
            end

        end

        function [obj]= checkIntensity(obj, r_missile_ECI, beta)
            % Function that checks if the visible intensity of the
            % satellite to an IR sensor is high enough to be detected
            % beta: current infrared output of missile based on plume (if
            % in boost), area observed, emissivity, etc.
            
            r_satellite_ECI = obj.x(1:3);
            r_sat_to_missile = r_missile_ECI - r_satellite_ECI;
            I = (obj.IR_sensor_gain*beta)/(norm(r_sat_to_missile)^2);

            if I >= obj.I_threshold
                obj.meetsIntensity = true;
            else
                obj.meetsIntensity = false;
            end

        end

        function [obj]= checkObservability(obj, r_missile_ECI, beta)
            % Overall function to check if a satellite can observe/detect
            % the missile, conditions are having line of sight and meeting
            % the intensity threshold
            obj.checkLOS(r_missile_ECI);
            obj.checkIntensity(r_missile_ECI, beta);

            if obj.hasLOS % && meetsIntensity
                obj.canObserve = true;
            else
                obj.canObserve = false;
            end

        end

        function [y]= getMeasurement(obj, x_missile_ECI, R)
            % Function to collect either a predicted measurement (R = 0)
            % or an actual measurement (R \neq 0)
            % Measurement vector: [az; el; az rate; el rate]

            x_satellite_ECI = obj.x;
            r_sat_to_missile = x_missile_ECI(1:3) - x_satellite_ECI(1:3);
            v_sat_to_missile = x_missile_ECI(4:6) - x_satellite_ECI(4:6);

            az = atan2(r_sat_to_missile(2), r_sat_to_missile(1)) + normrnd(0, sqrt(R(1,1)));
            el = atan2(-r_sat_to_missile(3), sqrt(r_sat_to_missile(1)^2 + r_sat_to_missile(2)^2)) + normrnd(0, sqrt(R(2,2)));
            

            az_rate = (r_sat_to_missile(1)*v_sat_to_missile(2) - r_sat_to_missile(2)*v_sat_to_missile(1))/(r_sat_to_missile(1)^2 + r_sat_to_missile(2)^2) + normrnd(0, sqrt(R(3,3)));

            el_rate = (r_sat_to_missile(3)*(r_sat_to_missile(1)*v_sat_to_missile(1) + r_sat_to_missile(2)*v_sat_to_missile(2)) - v_sat_to_missile(3)*(r_sat_to_missile(1)^2 + r_sat_to_missile(2)^2))/ ...
                (r_sat_to_missile(1)^2 + r_sat_to_missile(2)^2 + (r_sat_to_missile(3)^2)*sqrt(r_sat_to_missile(1)^2 + r_sat_to_missile(2)^2)) + normrnd(0, sqrt(R(4,4)));

            y = [az; el; az_rate; el_rate];

        end

        function [obj]= formHfunction(obj)
            % Function to form a matrix-valued function which, when
            % evaluated, forms the H matrix corresponding to the current
            % satellite's measurement

            syms x y z xdot ydot zdot real

            h1 = atan2(y, x);
            h2 = atan2(-z, sqrt(x^2 + y^2));         
            
            h3 = (x*ydot - y*xdot)/(x^2 + y^2);
            h4 = (z*(x*xdot + y*ydot) - zdot*(x^2 + y^2))/ ...
                (x^2 + y^2 + (z^2)*sqrt(x^2 + y^2));

            H_sym = jacobian([h1;h2;h3;h4], [x; y; z; xdot; ydot; zdot]);
            obj.H_func = matlabFunction(H_sym, "Vars", {x, y, z, xdot, ydot, zdot});
        end

        function [H]= evaluateH(obj, x_missile)
            obj.H = obj.H_func(x_missile(1), x_missile(2), x_missile(3), x_missile(4), x_missile(5), x_missile(6));
            H = obj.H;
        end
        
        function[r, v]= elements_to_rv(obj, a, e, i, O, w, f, mu)


            %Input all element angles in degrees
            %Input mu in SI units, no km stuff
            
            i = i*pi/180;
            O = O*pi/180;
            w = w*pi/180;
            f = f*pi/180;
            
            r_norm = (a*(1-e^2))/(1+e*cos(f));
            r_peri = [r_norm*cos(f); r_norm*sin(f); 0];
            
            p = a*(1-e^2);
            v_peri = [sqrt(mu/p)*-1*sin(f); sqrt(mu/p)*(e+cos(f)); 0];
            
            % Rotation matrix terms
            R_11 = cos(O)*cos(w) - sin(O)*sin(w)*cos(i);
            R_12 = -cos(O)*sin(w) - sin(O)*cos(w)*cos(i);
            R_13 = sin(O)*sin(i);
            
            R_21 = sin(O)*cos(w) + cos(O)*sin(w)*cos(i);
            R_22 = -sin(O)*sin(w) + cos(O)*cos(w)*cos(i);
            R_23 = -cos(O)*sin(i);
            
            R_31 = sin(w)*sin(i);
            R_32 = cos(w)*sin(i);
            R_33 = cos(i);
            
            R = [R_11, R_12, R_13;
                 R_21, R_22, R_23;
                 R_31, R_32, R_33];
            
            r = R*r_peri;
            v = R*v_peri;
            
         end




    end
end

