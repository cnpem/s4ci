import numpy as np

##############################################
##            LIPDIS STRUCTURES             ##
##############################################

##### Variables:
structures = ['Lamellar', 'Inverse Hexagonal', 'Primitive Cubic', 'Diamond Cubic', 'Gyroid Cubic']

##### Functions:

#### Classes:
class Structure:
    def __init__(self):
        #self.peaks = []
        return
    
class Lamellar(Structure):
    def __init__(self):
        super().__init__()
        self.peaks_spacings = np.array([1,2,3])
        self.name = 'Lamellar'

    def lattice_parameter(self, q, peak_spacings):
        self.a = 2*np.pi*peak_spacings[0]/q
        return self.a
    
    def peak_position(self, peak_spacings, a=False):
        #TODO ver melhor se a como passar somente algum argumento pelo nome e sem ordem especifica quando chamar a funcção
        if a:
            pass
        else:
            a = self.a

        q = 2*np.pi*peak_spacings/a
        return q
    
class InverseHexagonal(Structure):
    def __init__(self):
        super().__init__()
        self.peaks_spacings = np.array([1,np.sqrt(3),np.sqrt(4)])
        self.name = 'Inverse Hexagonal'

    def lattice_parameter(self, q, peak_spacings):
        self.a = (2/np.sqrt(3))*(2*np.pi/q)*peak_spacings[0]
        return self.a
    
    def peak_position(self, peak_spacings, a=False):
        if a:
            pass
        else:
            a = self.a

        q = (2/np.sqrt(3))*(2*np.pi/a)*peak_spacings
        return q
    
class Cubic(Structure):
    def __init__(self):
        super().__init__()
        #self.name = ''

    def lattice_parameter(self, q, peak_spacings):
        self.a = (2*np.pi/q)*peak_spacings[0]
        return self.a
        
    def peak_position(self, peak_spacings, a=False):
        if a:
            pass
        else:
            a = self.a
        q = (2*np.pi/a)*peak_spacings
        return q
        
class PrimitiveCubic(Cubic):
    def __init__(self):
        super().__init__()
        self.peaks_spacings = np.array([np.sqrt(2),np.sqrt(4),np.sqrt(6),np.sqrt(8),np.sqrt(10),np.sqrt(12),np.sqrt(14)])
        self.name = 'Primitive Cubic'

class DiamondCubic(Cubic):
    def __init__(self):
        super().__init__()
        self.peaks_spacings = np.array([np.sqrt(2),np.sqrt(3),np.sqrt(4),np.sqrt(6),np.sqrt(8),np.sqrt(9),np.sqrt(10), np.sqrt(11)])
        self.name = 'Diamond Cubic'
    
class GyroidCubic(Cubic):
    def __init__(self):
        super().__init__()
        self.peaks_spacings = np.array([np.sqrt(6),np.sqrt(8),np.sqrt(14),np.sqrt(16),np.sqrt(20),np.sqrt(22),np.sqrt(24)])
        self.name = 'Gyroid Cubic'

dic = {0:Lamellar, 1:InverseHexagonal, 2:PrimitiveCubic, 3:DiamondCubic, 4:GyroidCubic}
#todo  ver se esse arquivo e classes realmente fazem sentido
#todo  se sim continuar, senao pensar em nova estrategia para analise de picos e entao melhorar