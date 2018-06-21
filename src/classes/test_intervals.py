#!/usr/bin/env python

import cvxopt

class Test_ints:
    """
    Test Interval code
    """

    def __init__(self):
        # Define positive semi-definite matrix p
        P_low = cvxopt.matrix([8., -1., -1., 20.], (2, 2))
        a_low = cvxopt.matrix([-10., 2.], (2, 1))
        C_low = cvxopt.matrix([1., 3., -2., -4.], (2, 2))
        d_low = cvxopt.matrix([1., 4.], (2, 1))

        P_high = cvxopt.matrix([20., 1., 1., 40.], (2, 2))
        a_high = cvxopt.matrix([-6., 3.], (2, 1))
        C_high = cvxopt.matrix([2., 3., 8., 6.], (2, 2))
        d_high = cvxopt.matrix([10., 6.], (2, 1))

        cvxopt.solvers.options['show_progress'] = False
        sol = cvxopt.solvers.qp(P_low, a_low, C_low, d_high)
        print("Lower bounds: ")
        print(sol['x'])

        sol = cvxopt.solvers.qp(P_high, a_high, C_high, d_low)
        print("Upper bounds: ")
        print(sol['x'])