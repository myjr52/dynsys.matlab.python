clear;

% true probability density function (pdf)
var_x = 1;
mean_x = 0;
Omega_x = linspace(-5,5,1000);
px = (1/(sqrt(2*pi*var_x)))*exp(-(Omega_x-mean_x).^2/(2*var_x));

figure(1); clf;
plot(Omega_x,px,'LineWidth',2);
hold on;

% generate N random numbers with the mean zero and the variance 1 using
% randn
N_all = [100 10000];
x_bin = linspace(-5,5,30);
dx=mean(diff(x_bin));
line_style = {'rs-' 'go-'};
for idx=1:length(N_all)
    N_trial = N_all(idx);
    x_rand = randn(1,N_trial);
    
    % number of occurance of x_rand in x_bin
    N_occur = histcounts(x_rand,x_bin);
    
    figure(1);
    plot(x_bin(1:end-1)+dx/2, N_occur/(dx*N_trial),line_style{idx});
end

figure(1);
set(gca,'FontSize',14);
xlabel('Random Variable x Sampling Space: $\Omega_x$','Interpreter','latex');
ylabel('probability density function');
legend('True $p(x)$','N=100','N=10,000','Location','northeast','Interpreter','latex');