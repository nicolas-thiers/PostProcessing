# # Class definition

class mesh():
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []

class field():
    def __init__(self):
        pass

class volume():
    def __init__(self):
        self.history = str()
        self.time = []
        self.mesh = mesh()
        self.field = field()
