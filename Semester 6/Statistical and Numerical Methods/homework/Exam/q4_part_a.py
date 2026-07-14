import prettytable

from q4_common import (y1, y2, y3, load_excel, check_fit_quality,
                       fit_with_t_errors)

def part_a(file_path):
    worksheet_name = 'מדידות מייצגות'
    data = load_excel(file_path, worksheet_name)
    t = data['times']
    y = data['amplitudes']
    sigma_t = data['time_errors']
    sigma_y = data['amplitude_errors']
    results = {}
    for model, name in ((y1, "Y1"), (y2, "Y2"), (y3, "Y3")):
        fitted_params, chi_squared, sigma_eff = fit_with_t_errors(model, t, y, sigma_t, sigma_y)
        red, p = check_fit_quality(chi_squared, len(t), len(fitted_params))
        results[name] = (fitted_params, chi_squared, red, p)
        # plot_fit(model, t, y, sigma_y, sigma_t, fitted_params, f"Fit for {name}")
        # plot_chi2_landscape(model, t, y, sigma_y, name)
    table = prettytable.PrettyTable()
    table.field_names = ["Model", "A", "B", "C", "D", "Chi2", "Reduced Chi2", "p-value"]
    for name, (params, chi_squared, red, p) in results.items():
        table.add_row([name] + [f"{param:.4f}" for param in params] + [f"{chi_squared:.4f}", f"{red:.4f}", f"{p:.3g}"])
    print(table)

if __name__ == "__main__":
    part_a('Exam1.xlsx')
