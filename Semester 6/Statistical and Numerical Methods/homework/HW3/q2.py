import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def part1(pdf, x, domain=(-sp.oo, sp.oo)):

    moments = []
    for k in range(1, 5):
        moment = sp.integrate(x**k * pdf, (x, *domain))
        moments.append(sp.simplify(moment))
    return moments

def part2(dagum_pdf, moments, a, b, p, x_sym, b_vals, p_vals, a_values):
    # plot the pdf for different values of a
    x = np.linspace(0.01, 6, 100)
    plt.figure(figsize=(10, 6))
    for a_val in a_values:
        for b_val in b_vals:
            for p_val in p_vals:
                pdf_func = sp.lambdify(x_sym, dagum_pdf.subs({a: a_val, b: b_val, p: p_val}), 'numpy')
                plt.plot(x, pdf_func(x), label=f'a={a_val}, b={b_val}, p={p_val}')

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Dagum PDF for Different Values of a, b, and p')
    plt.legend()
    plt.grid(True)
    plt.show()

    moment_values = {}
    for a_val in a_values:
        for b_val in b_vals:
            for p_val in p_vals:
                for idx, moment in enumerate(moments):
                    identifier = f"a={a_val}, b={b_val}, p={p_val}, moment={idx + 1}"
                    moment = moments[idx].subs({b: b_val, p: p_val, a: a_val})
                    if sp.pi not in moment.free_symbols and not moment.has(sp.pi):
                        moment = moment.evalf(3)
                    else:
                        moment = sp.nsimplify(moment, [sp.pi], rational=False)
                        moment = sp.simplify(moment)
                    moment_values[identifier] = moment

    return moment_values

def part3(dagum_pdf, fourth_moment, a, b, p, x_sym):
    a_val = 10
    b_val = 0.01
    p_val = 1
    pdf_func = sp.lambdify(x_sym, dagum_pdf.subs({a: a_val, b: b_val, p: p_val}), 'numpy')
    
    print(f"Fourth moment for a={a_val}, b={b_val}, p={p_val}: {fourth_moment.subs({a: a_val, b: b_val, p: p_val})}")
    
    x = np.linspace(0.01, 6, 100)
    plt.figure(figsize=(10, 6))
    plt.plot(x, pdf_func(x), label=f'a={a_val}, b={b_val}, p={p_val}')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Dagum PDF for a=10, b=0.01, p=1')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    x = sp.Symbol('x', positive=True)
    p, a, b = sp.symbols('p a b', positive=True)

    dagum_pdf = (a * p / x) * (x / b)**(a*p) / (1 + (x / b)**a)**(p + 1)

    moments = part1(dagum_pdf, x, domain=(0, sp.oo))
    for i, m in enumerate(moments, 1):
        print(f"E[X^{i}] = {sp.latex(m)}")

    # b_vals, p_vals, a_vals = [1], [2], [0.5, 1, 2, 3, 4]
    # part_2_results = part2(dagum_pdf, moments, a, b, p, x, b_vals, p_vals, a_vals)
    # for identifier, moment in part_2_results.items():
    #     print(f"{identifier}: E[X] = {moment}")
    # print("\n")
    # b_vals, p_vals, a_vals = [0.5,1,2,3,4], [2], [4]
    # part_3_results = part2(dagum_pdf, moments, a, b, p, x, b_vals, p_vals, a_vals)
    # for identifier, moment in part_3_results.items():
    #     print(f"{identifier}: E[X] = {moment}")
    # print("\n")
    # b_vals, p_vals, a_vals = [1], [0.5,1,2,3,4], [4]
    # part_4_results = part2(dagum_pdf, moments, a, b, p, x, b_vals, p_vals, a_vals)
    # for identifier, moment in part_4_results.items():
    #     print(f"{identifier}: E[X] = {moment}")

    part3(dagum_pdf, moments[-1], a, b, p, x)

if __name__ == "__main__":
    main()