%2025/01/01
%LLG equation macrospin simulation
%Author: Zhang Tianyi
%Institute of Physics, Chinese Academy of Science
%E-mail: zhangty@iphy.ac.cn
clc
clear
%parameter 1.6mu_B = 1.5*10^-23 J T^-1
gamma = 3*power(10,11);%Hz T^-1
mu_s = 1.7*power(10,-22);%A m^2
K = 1.5*power(10,-23);%J
alpha = 0.01;
sigma_theta = pi/2;
sigma_phi = pi/2;
sigma = [sin(sigma_theta)*cos(sigma_phi),sin(sigma_theta)*sin(sigma_phi),cos(sigma_theta)];
variance_value = [];
variance_value_x = -0.36:0.02:0.36;
variance_value_x = variance_value_x';

var_n = [];
% n = 9.0;
% Hz = [-17.3:0.5:-15.3,-15.2:0.1:-13.0,-11.0:2:13.0,13.1:0.1:15.3];
for n = [-17.3:0.5:-15.3,-15.2:0.1:-13.0,-11.0:2:9.0]
    n
H_x = [];
mz_t = [];
my_t = [];
%ij=1.4not switch,1.5switch
% for ij = -2:0.05:2
for ij = n:0.005:n
%     ij
    mz = [];
    my = [];
    J_damp0 = 0*power(10,-23);
    J_field = 0;
    H_applied_value0 = ij*power(10,-24);
    H_x = [H_x,H_applied_value0];
    H_applied_theta = 0;
    H_applied_phi = 0;
    H_applied0 = H_applied_value0*[sin(H_applied_theta)*cos(H_applied_phi),sin(H_applied_theta)*sin(H_applied_phi),cos(H_applied_theta)];
    njj = 10;
    for jj  = 1:njj
%         jj
        %initial state
        m_x =[0.001];
        m_y =[0];
        m_z =[-sqrt(1-m_x(end)*m_x(end))];
        mx_end = m_x(end);
        my_end = m_y(end);
        mz_end = m_z(end);
        %time step and total time
        dt = 0.01*(mu_s*(1+alpha*alpha))/(abs(H_applied_value0)+abs(K)+abs(J_field)+abs(J_field))/gamma;
        N_t = 1000000;
        t = N_t*dt;
        progress = 1;
        w = 2*pi/njj/dt;
        %magnetic dynamics governed by llg equation
        for ii = 1:N_t
            %         if  10*ii/N_t >= progress
            %             progress = progress+1;
            %             fprintf('Progress is: %d%%\n', 100*ii/N_t);
            %         end
            J_damp = J_damp0;
            H_applied = H_applied0+0.03*power(10,-23)*sin(w*jj*dt)*[1,0,0];
            B_eff = [H_applied(1),H_applied(2),H_applied(3)+K*mz_end];
            S = [mx_end, my_end, mz_end];
            d = cross(S,B_eff+alpha*cross(S,B_eff))+J_damp*cross(S,cross(S,sigma)-alpha*sigma)+J_field*cross(S,sigma+alpha*cross(S,sigma));
            m_x_temp = mx_end-dt*gamma/(mu_s*(1+alpha*alpha))*d(1);
            m_y_temp = my_end-dt*gamma/(mu_s*(1+alpha*alpha))*d(2);
            m_z_temp = mz_end-dt*gamma/(mu_s*(1+alpha*alpha))*d(3);
            m_x_temp0 = m_x_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
            m_y_temp0 = m_y_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
            m_z_temp0 = m_z_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
            mx_end = m_x_temp0;
            my_end = m_y_temp0;
            mz_end = m_z_temp0;
%             if(ii > 0.8*N_t && mod(ii,100) == 1)
%                 m_x = [m_x,m_x_temp0];
%                 m_y = [m_y,m_y_temp0];
%                 m_z = [m_z,m_z_temp0];
%             end
%             if( mod(ii,100) == 1)
%                 m_x = [m_x,m_x_temp0];
%                 m_y = [m_y,m_y_temp0];
%                 m_z = [m_z,m_z_temp0];
%             end
        end

%         while var(m_z(end-200:end-100)) < var(m_z(end-100:end)) || abs(my_end) < 0.91
% %             abs(my_end)
%             for ii = 1:N_t
%                 J_damp = J_damp0;
%                 H_applied = 0.001*power(10,-22)*sin(w*jj*dt)*[sin(H_applied_theta)*cos(H_applied_phi),sin(H_applied_theta)*sin(H_applied_phi),cos(H_applied_theta)];
%                 B_eff = [H_applied(1),H_applied(2)+K*my_end,H_applied(3)];
%                 m_x_temp = mx_end-dt*gamma/(mu_s*(1+alpha*alpha))*(my_end*B_eff(3)-mz_end*B_eff(2)+alpha*((mx_end*B_eff(1)+my_end*B_eff(2)+mz_end*B_eff(3))*mx_end-B_eff(1))+J_damp*(mx_end*my_end*sigma(2)-my_end*my_end*sigma(1)-mz_end*mz_end*sigma(1)+mx_end*mz_end*sigma(3)-alpha*(-mz_end*sigma(2)+my_end*sigma(3)))+J_field*(-mz_end*sigma(2)+my_end*sigma(3)+alpha*(-my_end*my_end*sigma(1)-mz_end*mz_end*sigma(1)+mx_end*my_end*sigma(2)+mx_end*mz_end*sigma(3))));
%                 m_y_temp = my_end-dt*gamma/(mu_s*(1+alpha*alpha))*(mz_end*B_eff(1)-mx_end*B_eff(3)+alpha*((mx_end*B_eff(1)+my_end*B_eff(2)+mz_end*B_eff(3))*my_end-B_eff(2))+J_damp*(mx_end*my_end*sigma(1)-mx_end*mx_end*sigma(2)-mz_end*mz_end*sigma(2)+my_end*mz_end*sigma(3)-alpha*(-mx_end*sigma(3)+mz_end*sigma(1)))+J_field*(mz_end*sigma(1)-mx_end*sigma(3)+alpha*(mx_end*my_end*sigma(1)-mx_end*mx_end*sigma(2)-mz_end*mz_end*sigma(2)+my_end*mz_end*sigma(3))));
%                 m_z_temp = mz_end-dt*gamma/(mu_s*(1+alpha*alpha))*(mx_end*B_eff(2)-my_end*B_eff(1)+alpha*((mx_end*B_eff(1)+my_end*B_eff(2)+mz_end*B_eff(3))*mz_end-B_eff(3))+J_damp*(mx_end*mz_end*sigma(1)+my_end*mz_end*sigma(2)-mx_end*mx_end*sigma(3)-my_end*my_end*sigma(3)-alpha*(-my_end*sigma(1)+mx_end*sigma(2)))+J_field*(-my_end*sigma(1)+mx_end*sigma(2)+alpha*(mx_end*mz_end*sigma(1)+my_end*mz_end*sigma(2)-mx_end*mx_end*sigma(3)-my_end*my_end*sigma(3))));
%                 m_x_temp0 = m_x_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
%                 m_y_temp0 = m_y_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
%                 m_z_temp0 = m_z_temp/sqrt(m_x_temp*m_x_temp+m_y_temp*m_y_temp+m_z_temp*m_z_temp);
%                 mx_end = m_x_temp0;
%                 my_end = m_y_temp0;
%                 mz_end = m_z_temp0;
%                 if(ii > 0.8*N_t && mod(ii,100) == 1)
%                     m_x = [m_x,m_x_temp0];
%                     m_y = [m_y,m_y_temp0];
%                     m_z = [m_z,m_z_temp0];
%                 end
%             end
%         end
        mz = [mz, mz_end];
        my = [my, my_end];
    end
    variance_value = [variance_value, var(mz)];
    my_t = [my_t, my_end];
    mz_t = [mz_t, mz_end];
end
H_x0 = 0.03*power(10,-23);
Hz = [-17.3*power(10,-24):0.5*power(10,-24):-15.3*power(10,-24),-15.2*power(10,-24):0.1*power(10,-24):-13.0*power(10,-24),-11.0*power(10,-24):2*power(10,-24):13.0*power(10,-24),13.1*power(10,-24):0.1*power(10,-24):15.3*power(10,-24)];
var_ana = power(H_x0,4)./(32*power(Hz-K,4));
var_ana = var_ana';
% H_x = H_x';
% variance_value = variance_value';
% mz
mz_end
variance_value
var_n = [var_n,variance_value(end)];
end
var_n = var_n';
% plot(H_x,variance_value);
% my_t = my_t';
% plot(H_x,my_t);

% mx_end
% mz_end
% mz_end/mx_end
% K/J_damp

% x = 1:length(m_x);
% plot(x,m_x)
% hold on
% plot(x,m_y)
% hold on
% plot(x,m_z)
% ylim([-1,1])
% hold off
% legend("mx","my","mz")