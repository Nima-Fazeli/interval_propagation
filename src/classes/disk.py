#!/usr/bin/env python

import numpy as np
import cvxopt
import matplotlib.pyplot as plt
from matplotlib import animation

# import pdb
# animate globals
fig_animation = plt.figure(4)
ax_animation = plt.axes(xlim=(-2, 2), ylim=(-0.5, 2))
ax_animation.set_aspect('equal', adjustable='box')
ax_animation.grid()
line, = ax_animation.plot([], [], lw=2)
line.set_xdata([])
line.set_ydata([])


class Disk:
    """
    class for disk moving in the plane
    """

    def __init__(self):
        """
        Pre define the required variables for states and contact geometry
        """
        # simulation variables
        self.dt = 0.001
        self.simulate_time = 700
        self.g = np.array([0., -9.81, 0.])

        # geometry and inertia constants
        self.m = 1.
        self.R = 0.5
        self.rg = self.R * self.R / 2
        self.M = self.m * np.array([[1., 0., 0.], [0., 1., 0.], [0., 0.0, self.rg]])
        self.mu = 0.3
        self.ep = 0.5

        # initialize state variables
        self.q0 = np.array([0.0, 0.7, 0.0])
        self.v0 = np.array([0.1, 0.0, 0.0])
        self.q = self.q0
        self.v = self.v0
        self.v_next = np.array([0., 0., 0.0])
        self.J = np.array([[1., 0., self.R], [0., 1., 0.]])
        self.p = np.array([0., 0.])

        # initialize trajectories
        self.q_traj = self.q
        self.v_traj = self.v

        # operating intervals
        self.R_int = np.array([0.49, 0.51])
        self.m_int = np.array([0.99, 1.01])
        self.mu_int = np.array([0.20, 0.40])
        self.ep_int = np.array([0.45, 0.55])

        # figure handles
        plt.figure(1)
        self.fig0, self.ax0 = plt.subplots()
        plt.figure(2)
        self.fig1, self.ax1 = plt.subplots()
        plt.figure(3)
        self.fig2, self.ax2 = plt.subplots()
        # figure handles animation
        #self.fig_animation = plt.figure(4)
        #self.ax_animation = plt.axes(xlim=(-1, 1), ylim=(-2, 2))
        # self.line, = self.ax_animation.plot([], [], lw=2)

    def get_jacobian(self, bounds):
        if bounds == 'nominal':
            radius = self.R
        elif bounds == 'lower':
            radius = self.R_int[0]
        else:
            radius = self.R_int[1]
        jacobian = np.array([[1., 0., radius], [0., 1., 0.]])
        return jacobian

    @staticmethod
    def quadratic_optimization(a, b, a_c, b_c):
        # define linear constraints for qp
        b_const = np.vstack((b_c[:, None], np.zeros((4, 1))))
        a_const = np.vstack((-a_c, -np.eye(4)))

        # prepare variables for convex optimization
        p_opt = cvxopt.matrix(a + np.transpose(a))
        q_opt = cvxopt.matrix(b)
        a_const = cvxopt.matrix(a_const)
        b_const = cvxopt.matrix(b_const)

        # perform optimization
        cvxopt.solvers.options['show_progress'] = False
        sol = cvxopt.solvers.qp(p_opt, q_opt, a_const, b_const)

        p_jac = np.array([[-1., 1., 0., 0.], [0., 0., 1., 0.]])
        impulse = p_jac.dot(sol['x'])

        return impulse

    def resolve(self, bounds):
        if bounds == 'nominal':
            b_c = self.compute_b(bounds)
        elif bounds == 'lower':
            b_c = self.compute_b('upper')
        else:
            b_c = self.compute_b('lower')

        # compute the LCP matrix and vectors
        a = self.compute_a(bounds)
        b = self.compute_b(bounds)
        a_c = self.compute_a(bounds)
        self.p = self.quadratic_optimization(a, b, a_c, b_c)

    def dynamic_step(self, bounds):
        # setup the parameters depending on the bounds
        if bounds == 'nominal':
            contact_threshold = self.R - 1e-3
        elif bounds == 'lower':
            contact_threshold = self.R_int[0] - 1e-3
        else:
            contact_threshold = self.R_int[1] - 1e-3
        # get correct jacobian
        jacobian = self.get_jacobian(bounds)

        # update the world
        self.v_next = self.v + self.dt * self.g
        q_test = self.dt * self.v_next + self.q

        # if contact is detected, resolve it first
        if q_test[1] < contact_threshold:  # contact check
            self.resolve(bounds)
        else:
            self.p = np.array([[0.], [0.]])

        # re-update the world
        # here I play a trick to go from matrix to array (column)
        f_contact = np.transpose(np.transpose(jacobian).dot(self.p))[:][0]
        p_force = np.transpose(np.linalg.solve(self.M, f_contact))

        self.v = self.v + self.dt * self.g + p_force
        self.q = self.dt * self.v + self.q

        self.q_traj = np.vstack((self.q_traj, self.q))
        self.v_traj = np.vstack((self.v_traj, self.v))

    def simulate(self, simulate_list=['nominal', 'lower', 'upper'], with_figures=True, with_reset=True):
        # step through the simulation for nominal and the two bounds
        for bounds in range(len(simulate_list)):
            for step in range(self.simulate_time):
                self.dynamic_step(simulate_list[bounds])
            if with_figures:
                self.plot_x()
                self.plot_y()
                self.plot_t()
            if with_reset:
                self.reset_states()
        if with_figures:
            self.plot_figures()

    def animate_loop(self):
        # bounds = 'nominal'
        self.simulate(simulate_list=['nominal'], with_figures=False, with_reset=False)
        q_nom = self.q_traj
        self.reset_states()
        self.simulate(simulate_list=['lower'], with_figures=False, with_reset=False)
        q_low = self.q_traj

        anim = animation.FuncAnimation(fig_animation, self.animate, fargs=(q_nom, q_low), frames=500,
                                       interval=20, blit=True)
        anim.save('sample.mp4', writer='ffmpeg')
        # plt.show()

    def animate(self, i, q_nom, q_low):
        # step simulation
        self.dynamic_step('nominal')
        # draw new state
        theta = np.linspace(0, 2*np.pi, 1000)
        x = q_nom[i][0] + self.R * np.cos(theta)
        y = q_nom[i][1] + self.R * np.sin(theta)
        line.set_data(x, y)
        return line,

    def compute_b(self, bounds):
        if bounds == 'nominal':
            restitution_coefficient = self.ep
            jacobian = self.get_jacobian(bounds)
            J_bar = np.vstack((-jacobian[0][:], jacobian))
        elif bounds == 'lower':
            restitution_coefficient = self.ep_int[1]
            jacobian_bottom = self.get_jacobian(bounds)
            jacobian_top = self.get_jacobian('upper')
            J_bar = np.vstack((-jacobian_top[0][:], jacobian_bottom))
        else:
            restitution_coefficient = self.ep_int[0]
            jacobian_bottom = self.get_jacobian(bounds)
            jacobian_top = self.get_jacobian('lower')
            J_bar = np.vstack((-jacobian_top[0][:], jacobian_bottom))

        # predefine the affine part of the LCP
        b1 = J_bar[0, :].dot(self.v + self.dt * self.g)
        b2 = J_bar[1, :].dot(self.v + self.dt * self.g)
        bn = J_bar[2, :].dot(self.v * (1 + restitution_coefficient) + self.dt * self.g)
        b = np.array([b1, b2, bn, 0.])

        return b

    def compute_a(self, bounds):
        # predefine the necessary components of the A matrix for LCP
        if bounds == 'nominal':
            friction_coefficient = self.mu
            jacobian = self.get_jacobian(bounds)
            J_bar = np.vstack((-jacobian[0][:], jacobian))
        elif bounds == 'lower':
            jacobian_bottom = self.get_jacobian(bounds)
            jacobian_top = self.get_jacobian('upper')
            J_bar = np.vstack((-jacobian_top[0][:], jacobian_bottom))
            friction_coefficient = self.mu_int[0]
        else:
            jacobian_bottom = self.get_jacobian(bounds)
            jacobian_top = self.get_jacobian('lower')
            J_bar = np.vstack((-jacobian_top[0][:], jacobian_bottom))
            friction_coefficient = self.mu_int[1]

        contact_mass = J_bar.dot(np.linalg.solve(self.M, np.transpose(J_bar)))
        e_vert = np.array([1., 1., 0., 0.])
        e_horz = np.array([-1., -1., friction_coefficient])

        # format the matrix correctly
        a_top = np.array(np.vstack((contact_mass, np.transpose(e_horz))))
        a = np.hstack((a_top, e_vert[:, None]))

        return a

    def reset_states(self):
        self.q = self.q0
        self.v = self.v0
        self.q_traj = self.q0
        self.v_traj = self.v0

    def plot_x(self):
        #plt.figure(1)
        self.ax0.plot(self.q_traj[:, 0])
        self.ax0.set(xlabel='increment', ylabel='x (m)', title='disk x test')
        self.ax0.grid()

    def plot_y(self):
        #plt.figure(2)
        self.ax1.plot(self.q_traj[:, 1])
        self.ax1.set(xlabel='increment', ylabel='y (m)', title='disk y test')
        self.ax1.grid()

    def plot_t(self):
        #plt.figure(3)
        self.ax2.plot(self.q_traj[:, 2])
        self.ax2.set(xlabel='increment', ylabel='theta (radians)', title='disk theta test')
        self.ax2.grid()

    def plot_figures(self):
        self.fig0.legend(('Nominal Parameters', 'Lower Bound Parameters', 'Upper Bound Parameters'), loc='upper right')
        self.fig1.legend(('Nominal Parameters', 'Lower Bound Parameters', 'Upper Bound Parameters'), loc='upper right')
        self.fig2.legend(('Nominal Parameters', 'Lower Bound Parameters', 'Upper Bound Parameters'), loc='upper right')
        plt.show()
