% flag_perf_tests

ell = [ 4 8 16 32 64 128 256 512];
accuracy = [ ...
    1.15e-15 ... % 4
    9.31e-15 ... % 8
    3.86e-14 ... % 16
    4.0e-13 ... % 32
    1.10e-12 ... % 64
    2.8e-12 ...  % 128
    2.0e-11 ... % 256
    2.0e-10 % 512
    ];
speed = [ ...
    0.00025 ...  % 4
    0.0016 ... % 8
    0.01 ... % 16
    0.07 ... % 32
    0.9 ... % 64
    12 ... % 128
    190 ... % 256
    2750 % 512
    ];

nb_samples = 2*ell.^3;
nb_samples_pow = str2mat('t04', 't07', 't10', 't13', 't16', 't19', 't22', 't25', 't28');

figure('Position',[100 100 700 700])

subplot(2,1,1)
loglog( nb_samples, (1e-14)*ell.^2, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, accuracy, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x1','FontSize',20)
ylabel('y1','FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
set(gca,'YTick',[1e-16 1e-14 1e-12  1e-10 1e-8] ,'FontSize',20)
axis([80 1.5*nb_samples(8) 1e-16 1e-8])

%title('Accuracy of overall transform','FontSize',20)


subplot(2,1,2)
loglog( nb_samples, (5e-6)*ell.^4, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, speed, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x2','FontSize',20)
ylabel('y2','FontSize',20)
axis([80 1.5*nb_samples(8) 1e-5 1e6])
set(gca,'YTick',[1e-6 1e-4 1e-2  1e0 1e2 1e4 1e6] ,'FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
%title('Speed of overall transform','FontSize',20)

