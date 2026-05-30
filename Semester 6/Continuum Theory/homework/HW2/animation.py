import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def compute_trajectory(starting_point, velocity_field, time_steps):
    trajectory = [starting_point]
    for t in time_steps[1:]:
        current_position = trajectory[-1]
        velocity = velocity_field(current_position, t)
        new_position = current_position + velocity * (time_steps[1] - time_steps[0])
        trajectory.append(new_position)
    return np.array(trajectory)

def animate_trajectory(trajectory, time_steps):
    fig, ax = plt.subplots()
    ax.set_xlim(-1, 10)
    ax.set_ylim(-1, 10)

    def update(frame):
        ax.cla()
        ax.set_xlim(-1, 10)
        ax.set_ylim(-1, 10)
        t = time_steps[frame]
        X, Y, U, V = prepare_streamlines_for_animation(
            lambda pos, _t=t: np.array([pos[1] * _t, 1]),
            (-1, 10), (-1, 10)
        )
        ax.streamplot(X, Y, U, V, color='blue', density=1.5)
        ax.plot(trajectory[:frame + 1, 0], trajectory[:frame + 1, 1], 'ro-')
        return []

    ani = animation.FuncAnimation(fig, update, frames=len(time_steps), blit=False)
    ani.save("trajectory_animation.mp4", writer='ffmpeg', fps=5)


def prepare_streamlines_for_animation(velocity_field, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 20)
    y = np.linspace(y_range[0], y_range[1], 20)
    X, Y = np.meshgrid(x, y)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            pos = np.array([X[i, j], Y[i, j]])
            velocity = velocity_field(pos)
            U[i, j] = velocity[0]
            V[i, j] = velocity[1]

    return X, Y, U, V

def main():
    starting_point = np.array([0, 0])
    alpha, beta = 1, 1
    velocity_field = lambda pos, t: np.array([alpha * pos[1] * t, beta])
    time_steps = np.linspace(0, 3, 15)
    trajectory = compute_trajectory(starting_point, velocity_field, time_steps)
    animate_trajectory(trajectory, time_steps)

if __name__ == "__main__":
    main()