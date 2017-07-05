from calculators.Calculator import Calculator

class AnaPot(Calculator):

    def __init__(self): 
        super(AnaPot, self).__init__()

    def get_energy(self, x, y):
        return 4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y

    def get_grad(self, x, y):
        return (
            4.5 + 2*x -2*y + 4*x**3 - 4*x*y,
            -4 + 4*y - 2*x - 2*x**2
        )

    def get_hessian(self, x, y):
        return ((12*x**2 + 2 - 4*y, -4*x-2),
                (-4*x-2, 4)
        )
