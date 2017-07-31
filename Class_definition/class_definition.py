# # Class definition

class mesh():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

class field():
    def __init__(self):
        self.U = []
        self.V = []
        self.W  = []
        self.P  = []
        self.T  = []
        
class gradient():
    def __init__(self):
        self.dUdx = []
        self.dUdy = []
        self.dUdz = []
        self.dVdx = []
        self.dVdy = []
        self.dVdz = []
        self.dWdx = []
        self.dWdy = []
        self.dWdz = []
        self.dPdx = []
        self.dPdy = []
        self.dPdz = []
        self.dTdx = []
        self.dTdy = []
        self.dTdz = []
        
class volume():
    def __init__(self):
        self.history = str()
        self.time = []
        self.mesh = mesh()
        self.field = field()
        self.gradient = gradient()
