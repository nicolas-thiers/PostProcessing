# # Class definition

class mesh():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

class field():
    def __init__(self):
        '''
		self.U = []
        self.V = []
        self.W  = []
        self.P  = []
        self.T  = []
        self.uu = []
        self.uv = []
        self.uw = []
        self.vv = []
        self.vw = []
        self.ww = []
        self.PP = []
        self.TT = []
        self.uT = []
        self.vT = []
        self.wT = []
        self.k  = []
        self.uuu  = []
        self.vvv  = []
        self.www  = []
        self.TTT  = []
        self.uuuu = []
        self.vvvv = []
        self.wwww = []
        self.TTTT = []
        self.I2 = []
        self.I3 = []
        self.QA = []
        self.RA = []
        self.QS = []
        self.RS = []
        self.QW = []
        self.RW = []
        '''
class gradient():
    def __init__(self):
        '''
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
        '''
class volume():
    def __init__(self):
        self.history = str()
        self.time = []
        self.mesh = mesh()
        self.field = field()
        self.gradient = gradient()
