clear all;
rng(4);

%% Parameters
n = 2; % Number of inputs/outputs
p = 1; % Type of norm for input constraint
a = 0.1; % Constraint on input norm

%% System
Q = rand(n);
for i = 1:n
    Q(:,i) = Q(:,i)/norm(Q(:,i));
end
D = diag(-sort(rand(n,1),'descend'));
A = Q*D*inv(Q);
B = rand(n);

f = @(y) A*y;
g = @(y) B;

dy = @(y,u) f(y) + g(y)*u;
hi = @(y) -inv(g(y))*f(y);

%% Control law
u_ref = 0*ones(n,1);
alpha = @(y) (a/norm(hi(y),p)); % This only works for u_ref=0
u = @(y) (1+alpha(y))*u_ref-alpha(y)*hi(y);

%% Simulation 
opt    = odeset('Events', @isConverged);

y0 = rand(n,1);

t = linspace(0,20,10000);
[t_step,y_step] = ode45(@(t,y) dy(y,u_ref), t, y0, opt);
[t_cont,y_cont] = ode45(@(t,y) dy(y,u(y)), t, y0, opt);

%% Compute the input
u_cont = zeros(size(y_cont));
for i = 1:length(t_cont)
    u_cont(i,:) = u(y_cont(i,:)')';
end

% Compute the constraint region
[U1,U2] = meshgrid(linspace(-2*a,2*a,200));
for i = 1:size(U1,1)
    for j = 1:size(U1,2)
        if norm([U1(i,j), U2(i,j)],p) > a
            U1(i,j) = 0;
            U2(i,j) = 0;
        end
    end
end

k = convhull([U1(:),U2(:)]);

%% Figures
figure(1);
clf; 
subplot(2,2,[1,2])
grid on; hold on;

plot([t_step; t(end)],[vecnorm(y_step',p), 0],'b','linewidth',2, 'DisplayName', '$y_{step}$')
plot([t_cont; t(end)],[vecnorm(y_cont',p), 0],'r--','linewidth',2, 'DisplayName', '$y_{cont}$')

title('Norm of step response','interpreter','latex', 'fontsize',18)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$\|y(t)\|$','interpreter','latex','fontsize',18)
legend('Interpreter','latex','fontsize',14)

subplot(2,2,3)
grid on; hold on;
plot(y_step(:,1),y_step(:,2),'b','linewidth',2, 'DisplayName', '$y_{step}$')
plot(y_cont(:,1),y_cont(:,2),'r--','linewidth',2, 'DisplayName', '$y_{cont}$')
plot(y0(1), y0(2), 'gx', 'linewidth', 2, 'markersize', 20, 'DisplayName', '$y_{0}$')
plot(0,0, 'kx', 'linewidth', 2, 'markersize', 20, 'DisplayName', '$y_{ref}$')

title('Phase plot for output','interpreter','latex', 'fontsize',18)
xlabel('$y_1$','interpreter','latex','fontsize',18)
ylabel('$y_2$','interpreter','latex','fontsize',18)
legend('Interpreter','latex','fontsize',14,'location','northwest')

xlim([-1,1])
ylim([-1,1])

subplot(2,2,4)
grid on; hold on; daspect([1,1,1])
xlim([-1.5*a,1.5*a])
ylim([-1.5*a,1.5*a])
plot(u_cont(:,1),u_cont(:,2),'r--','linewidth',2,'DisplayName','$u_{cont}$')
fill(U1(k),U2(k), 'r','FaceAlpha',0.3, 'edgecolor','none','DisplayName','$\mathcal{U}$')
plot(0,0, 'kx', 'linewidth', 2, 'markersize', 20,'DisplayName','$u_{ref}$')

title('Phase plot for input','interpreter','latex', 'fontsize',18)
xlabel('$u_1$','interpreter','latex','fontsize',18)
ylabel('$u_2$','interpreter','latex','fontsize',18)
legend('Interpreter','latex','fontsize',14)

%% Functions
function [value, isterminal, direction] = isConverged(t, y)
    value      = (norm(y) < 1e-3);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end