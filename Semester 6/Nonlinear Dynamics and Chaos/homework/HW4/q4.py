import numpy as np
import matplotlib.pyplot as plt
import scienceplots as sp
plt.style.use(['science', 'notebook', 'grid'])

def main():
    K = 0.5
    delta = K**2 / 4  # bifurcation edge of the 1:2 tongue
    eta = np.pi + delta
    G = lambda theta: theta + eta - K * np.sin(theta)

    thetas = np.linspace(0, 2*np.pi, 1000)
    G2 = np.mod(G(G(thetas)), 2*np.pi)  # exact second iterate

    plt.plot(thetas, G2, label=r'$G^{(2)}(\theta) mod 2\pi$')
    plt.plot(thetas, thetas, label=r'diagonal $(\theta,\theta)$', color='red')
    plt.axvline(np.pi/4, linestyle='--', color='gray', label=r'predicted $\theta^*=\pi/4,\,5\pi/4$')
    plt.axvline(5*np.pi/4, linestyle='--', color='gray')
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\theta_2$')
    plt.title(f'Second iterate of the sine-circle map, K={K}, $\\eta=\\pi+K^2/4$')
    plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 3*np.pi/2, 7*np.pi/4, 2*np.pi],
               [r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$', r'$2\pi$'])
    plt.yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 3*np.pi/2, 7*np.pi/4, 2*np.pi],
               [r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$', r'$2\pi$'])
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
