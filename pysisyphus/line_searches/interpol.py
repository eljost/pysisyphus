import numpy as np

def interpol_alpha_quad(f_0, df_0, f_alpha_0, alpha_0):
    return -df_0*alpha_0**2 / 2 / (f_alpha_0 - f_0 - df_0*alpha_0)


def interpol_alpha_cubic(f_0, df_0, f_alpha_0, f_alpha_1, alpha_0, alpha_1):
    quot = 1 / (alpha_0**2 * alpha_1**2 * (alpha_1 - alpha_0))
    A = np.array(((alpha_0**2, -alpha_1**2),
                  (-alpha_0**3, alpha_1**3)))
    B = np.array(((f_alpha_1 - f_0 - df_0*alpha_1),
                  (f_alpha_0 - f_0 - df_0*alpha_0)))
    a, b = quot * A @ B
    alpha_cubic = (-b + (b**2 - 3*a*df_0)**(0.5)) / (3*a)
    return alpha_cubic
