function plot_comparison_res(info)
%PLOT_COMPARISON_RES Summary of this function goes here
%   Detailed explanation goes here
figure();

fig = figure();
t = tiledlayout(1, 2,'TileSpacing','Compact');

nexttile
x = [info.iter];
semilogy(x, [info.res_norm], 'DisplayName', 'QR fact');
hold on
semilogy(x, [info.res_norm_hutch1], 'DisplayName', 'Hutch++ 1');
semilogy(x, [info.res_norm_hutch5], 'DisplayName', 'Hutch++ 5');
semilogy(x, [info.res_norm_hutch10], 'DisplayName', 'Hutch++ 10');
semilogy(x, [info.res_norm_hutch25], 'DisplayName', 'Hutch++ 25');
semilogy(x, [info.res_norm_eigs], 'DisplayName', 'eigs');
semilogy(x, [info.res_norm_eigs5], 'DisplayName', 'eigs 5', 'Color', 'black');
semilogy(x, [info.res_norm_eigs25], 'DisplayName', 'eigs 25', 'Color', 'green');
ylabel("Residual")

nexttile
semilogy(x, [info.time_res_norm], 'DisplayName', 'QR fact');
hold on
semilogy(x, [info.res_norm_hutch1_time], 'DisplayName', 'Hutch++ 1');
semilogy(x, [info.res_norm_hutch5_time], 'DisplayName', 'Hutch++ 5');
semilogy(x, [info.res_norm_hutch10_time], 'DisplayName', 'Hutch++ 10');
semilogy(x, [info.res_norm_hutch25_time], 'DisplayName', 'Hutch++ 25');
semilogy(x, [info.res_norm_eigs_time], 'DisplayName', 'eigs');
semilogy(x, [info.res_norm_eigs5_time], 'DisplayName', 'eigs 5', 'Color', 'black');
semilogy(x, [info.res_norm_eigs25_time], 'DisplayName', 'eigs 25', 'Color', 'green');
semilogy(x, [info.time_assembly_res], 'DisplayName', 'Assembly', 'Color', 'blue')
legend()
ylabel("Time")

end

