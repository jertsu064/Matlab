function [thetas,vels,accs] = slidercrank_soln(r1,r2,r3,r4,theta2, theta2dot, theta2ddot, init_values)

% there are three nested functions - one each for position, velocity and
% acceleration analysis.
% Solve for theta3 and theta4
thetas = fsolve(@position_analysis,init_values);
% Use the angles from above for the velocity analysis
vels = velocity_analysis(thetas);
% Use the angles and velocities for the acceleration analysis
accs = acceleration_analysis(thetas,vels);
% projs=vec_proj(vector,cood);

    function F = position_analysis(v)
         theta3 = v(1);
         theta4 = v(2);
         F = [-r1+r2*cos(theta2)+r3*cos(theta3)-r4*cos(theta4);
             r2*sin(theta2)+r3*sin(theta3)-r4*sin(theta4)];
    end

    function F = velocity_analysis(thetas)
        theta3 = thetas(1);
        theta4 = thetas(2);
        A = [-r3*sin(theta3) r4*sin(theta4);r3*cos(theta3) -r4*cos(theta4)];
        B = [r2*sin(theta2)*theta2dot; -r2*cos(theta2)*theta2dot];
        F = A\B;
    end

    function F = acceleration_analysis(thetas,vels)
        theta3 = thetas(1);
        theta4=thetas(2);
%         r1dot = vels(1);
        theta3dot = vels(1);
        theta4dot = vels(2);
        A=[-r4*sin(theta4) -r4*sin(theta3) r3*sin(theta3) r3*cos(theta3); r4*cos(theta3) -r4*sin(theta4) -r3*cos(theta3) r3*sin(theta3)];
        B=[-r2*theta2ddot*sin(theta2) -r2*theta2dot.^2*cos(theta2);r2*theta2ddot*cos(theta2) -r2*theta2dot.^2*sin(theta4)];
        F=A\B;
    end

end
