reduction_vals = cat(2, zeros(1), linspace(0,100,1000));
wide_density = betapdf(0.01*reduction_vals,4.5,1.5);
normal_density = betapdf(0.01*reduction_vals,12,3);

figure(1)
clf(1)
plot(reduction_vals,wide_density,'DisplayName', 'Wide prior','LineWidth',3)
hold on
plot(reduction_vals, normal_density,'DisplayName','Narrow prior','LineWidth',3)
hold off
legend('FontSize',14)
xlabel('Percentage Reduction','FontSize',14)
ylabel('Probability density','FontSize',14)


[p_td_vals_1, ptd_pdf_normal] = Get_PDF_ptd(reduction_vals, normal_density,PreCalc);
[p_td_vals_2, ptd_pdf_wide] = Get_PDF_ptd(reduction_vals, wide_density,PreCalc);

figure(2)
clf(2)
hold on
plot(p_td_vals_2, ptd_pdf_wide,'DisplayName','Wide prior','LineWidth',3)
plot(p_td_vals_1,ptd_pdf_normal,'DisplayName','Narrow prior','LineWidth',3)
hold off
legend('FontSize',14)
xlabel('p target die','FontSize',14)
ylabel('Probability density','FontSize',14)