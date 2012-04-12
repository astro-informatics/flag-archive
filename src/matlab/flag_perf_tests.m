% flag_perf_tests

ell = [ 4 8 16 32 64 96 128 192 256];
accuracy = [ ...
    1.35e-15 ... % 4
    1.10e-14 ... % 8
    3.70e-14 ... % 16
    4.26e-13 ... % 32
    1.10e-12 ... % 64
    1.75e-12 ... % 96
    2.9e-12 ...  % 128
    2.0e-11 ...  % 192 
    3.44e-11 ... % 256
    ];
speed = [ ...
    0.0002 ...  % 4
    0.002 ... % 8
    0.01 ... % 16
    0.1 ... % 32
    1.2 ... % 64
    5 ...  % 96
    15 ... % 128
    70 ... % 192
    220 % 256
    ];

figure('Position',[100 100 600 600])

subplot(2,1,1)
loglog( ell, (1e-15)*ell.^2, 'red','LineWidth', 2 )
hold on
loglog(ell, accuracy, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
set(gca,'XTick',[  16 64  128  256 512 ],'FontSize',20)
set(gca,'YTick',[1e-16 1e-14 1e-12  1e-10 1e-8] ,'FontSize',20)
axis([0 512 1e-16 1e-8])
xlabel('N = L','FontSize',20)
ylabel('Average max error','FontSize',20)
title('Accuracy of overall transform','FontSize',20)


subplot(2,1,2)
loglog( ell, (1e-6)*ell.^4, 'red','LineWidth', 2 )
hold on
loglog(ell, speed, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
axis([0 512 1e-5 1e5])
xlabel('N = L','FontSize',20)
set(gca,'YTick',[1e-6 1e-4 1e-2  1e0 1e2 1e4 1e6] ,'FontSize',20)
ylabel('Average speed (s)','FontSize',20)
set(gca,'XTick',[  16 64  128  256 512 ],'FontSize',20)
title('Speed of overall transform','FontSize',20)

