import numpy as np
import math

class Solution(object):
    def __init__(self, height, tao_0, solid, fluid, gravity):
        self.height = height;
        self.tao_0 = tao_0;
        self.g = gravity;

        fluid_name = list(fluid.keys())[0]
        self.mu = fluid[fluid_name]["Viscosity"]
        self.c_f = fluid[fluid_name]["Compressibility"]
        self.rho_f = fluid[fluid_name]["Density"]

        solid_name = list(solid.keys())[0]
        self.c_s = solid[solid_name]["Compressibility"]
        self.permeability = solid[solid_name]["Permeability"]
        self.phi = solid[solid_name]["Porosity"]
        self.nu = solid[solid_name]["PoissonsRatio"]
        self.G = solid[solid_name]["ShearModulus"]
        self.rho_s = solid[solid_name]["Density"]


        self.rho = self.rho_f * self.phi + (1 - self.phi) * self.rho_s;
        self.lame1st = 2 * self.G * self.nu / ( 1 - 2 * self.nu );
        self.M = 2 * self.G + self.lame1st;

        self._calculate_K_s()
        self._calculate_K_f()
        self._calculate_K_phi()
        self._calculate_K()

        self._calculate_alpha()
        self._calculate_Q()

        self._calculate_H()
        self._calculate_R()
        self._calculate_K_p()
        self._calculate_B()
        self._calculate_K_u()
        self._calculate_ni_u()
        self._calculate_K_ni_u()
        self._calculate_gama()
        self._calculate_c()

        self._calculate_K_v_u()
        self._calculate_c_m()

    # Internal functions --------------------------------------------------------------------

    def _calculate_K_s( self ):
        if self.c_s == 0.0:
            self.K_s = 1.0e+100
        else:
            self.K_s = 1.0 / self.c_s

    def _calculate_K_f( self ):
        if self.c_f == 0.0:
            self.K_f = 1.0e+100
        else:
            self.K_f = 1.0 / self.c_f

    def _calculate_K_phi( self ):
        self.K_phi = self.K_s

    def _calculate_K( self ):
        self.K = 2*self.G*( 1 + self.nu ) / ( 3 - 6*self.nu )

    def _calculate_H( self ):
        self.H = 1.0 / ( ( 1.0 / self.K ) - ( 1.0 / self.K_s ) )

    def _calculate_R( self ):
        self.R = 1.0 / ( ( 1.0 / self.H ) + self.phi * ( ( 1.0 / self.K_f ) - ( 1.0 / self.K_phi ) ) )

    def _calculate_K_p( self ):
        self.K_p = self.phi * self.K / self.alpha

    def _calculate_B( self ):
        self.B = self.R / self.H

    def _calculate_K_u( self ):
        if self.alpha * self.B != 1.0:
            self.K_u = self.K / ( 1.0 - self.alpha * self.B )
        else:
            self.K_u = 1.0e+100

    def _calculate_ni_u( self ):
        self.nu_u = ( ( 3.0 * self.nu + self.alpha * self.B * ( 1.0 - 2.0 * self.nu ) ) /
                          ( 3.0 - self.alpha * self.B * ( 1.0 - 2.0 * self.nu ) ) )

    def _calculate_K_ni_u( self ):
        self.K_ni_u = ( 3.0 * self.K_u * ( 1.0 - self.nu_u ) ) / ( 1.0 + self.nu_u )

    def _calculate_gama( self ):
        self.gama = ( self.B * ( 1.0 + self.nu_u ) ) / ( 3.0 * ( 1.0 - self.nu_u ) )

    def _calculate_c( self ):
        self.c = ( ( 2.0 * self.permeability * self.G * ( 1.0 - self.nu ) * ( self.nu_u - self.nu ) ) /
                    ( self.mu * ( self.alpha ** 2.0 ) * ( 1.0 - self.nu_u ) * ( ( 1.0 - 2.0 * self.nu ) ** 2.0 ) ) )

    def _calculate_K_v_u( self ):
        self.K_v_u = 3*self.K_u*( 1 - self.nu_u ) / ( 1 + self.nu_u )

    def _calculate_c_m( self ):
        self.c_m = ( self.alpha * ( 1.0 + self.nu ) ) / ( 3.0 * self.K * ( 1.0 - self.nu ) )

    def _calculate_ni( self ):
        self.nu = ( 3.0 * self.K - 2.0 * self.G ) / ( 2.0 * ( 3.0 * self.K + self.G ) )

    def _calculate_alpha( self ):
        self.alpha = 1.0 - ( self.K / self.K_s )

    def _calculate_Q( self ):
        self.Q = 1 / ( self.phi * self.c_f + ( self.alpha - self.phi ) * self.c_s )

    def _calculate_p_0( self, yPosition ):
        self.p_0 = self.alpha * self.Q / ( self.M + self.alpha * self.alpha * self.Q ) * \
            ( self.tao_0 + 0.5 * self.rho * self.g * self.height ) - self.rho_f * self.g * \
            ( yPosition - 0.5 * self.height );
        return self.p_0

    def _calculate_v_0( self, yPosition ):
        self.v_0 = ( self.rho - self.alpha * self.rho_f ) * self.g / ( 2 * self.M ) * yPosition \
            * yPosition - ( ( self.tao_0 + 0.5 * self.rho * self.g * self.height ) / ( self.M + \
            self.alpha * self.alpha * self.Q ) + ( self.rho - self.alpha * self.rho_f ) * self.g \
            *self.height / ( 2 * self.M ) ) * yPosition;
        return self.v_0

    def __getFunctionValuesAtTime(self, function, time, ny, numberOfSummationTerms):
        positionValues = self.getPositionValues(ny)
        functionValues = np.zeros(len(positionValues))
        for i,pos in enumerate(positionValues):
            functionValues[i] = function(pos, time, numberOfSummationTerms)
        return functionValues

    def __getFunctionValuesAtPosition(self, function, position, timeList, numberOfSummationTerms):
        functionValues = np.zeros(len(timeList))
        for i,t in enumerate(timeList):
            functionValues[i] = function(position, t, numberOfSummationTerms)
        return functionValues




    # Class interface -----------------------------------------------------------------------
    def getPositionValues(self, ny=200):
        '''If ny is an interger, then a list of equally spaced vertical position will be created. However, ny can be
            a list or an numpy array with specified vertical positions.'''
        if type(ny) == int:
            return np.linspace(0, self.height, ny)
        elif type(ny)==np.ndarray or type(ny) == list:
            return ny
        else:
            raise

    def getPressureValuesConstTime(self, time, ny=200, numberOfSummationTerms=200):
        return self.__getFunctionValuesAtTime(self.getPressureValue, time, ny, numberOfSummationTerms)

    def getDisplacementValuesConstTime(self, time, ny=200, numberOfSummationTerms=200):
        return self.__getFunctionValuesAtTime(self.getDisplacementValue, time, ny, numberOfSummationTerms)

    def getPressureValuesAtPosition(self, position, timeList, numberOfSummationTerms=200):
        return self.__getFunctionValuesAtPosition(self.getPressureValue, position, timeList, numberOfSummationTerms)

    def getDisplacementValuesAtPosition(self, position, timeList, numberOfSummationTerms=200):
        return self.__getFunctionValuesAtPosition(self.getDisplacementValue, position, timeList, numberOfSummationTerms)


    def getPressureValue(self, yPosition, time, numberOfSummationTerms=200):
        position = yPosition;
        if time == 0.0:
            pressureValue = self._calculate_p_0( yPosition )
            return pressureValue
        else:
            summationResult = 0
            for j in range(1, numberOfSummationTerms):
                term_1 = (-1.0)**(j+1)/(2*j - 1)
                term_2 = math.exp( -(2*j - 1)**2.0 * time*self.c*math.pi**2.0 / (4.0*self.height**2.0) )
                term_3 = math.cos( (2*j - 1)*math.pi*position / (2.0*self.height) )
                summationResult += term_1*term_2*term_3
            barP0 = self.alpha * self.Q / ( self.M + self.alpha * self.alpha * self.Q ) * \
                ( self.tao_0 + 0.5 * self.rho * self.g * self.height ) - self.rho_f * self.g * \
                ( 0.5 * self.height );
            pressureValue = 4.0 * barP0 * summationResult / math.pi
            pressureValue += self.rho_f * self.g * ( self.height - yPosition )
            return pressureValue


    def getDisplacementValue(self, yPosition, time, numberOfSummationTerms=200):
        # position = self.height - yPosition;
        position = yPosition;
        if time == 0.0:
            displacementValue = self._calculate_v_0( yPosition )
            return displacementValue
        else:
            summationResult = 0
            for j in range(1, numberOfSummationTerms):
                term_1 = (-1)**(j+1) / (2*j - 1)**2
                # term_1 = 1 / (2*j - 1)**2
                term_2 = math.exp( -(2*j - 1)**2 * self.c*time*math.pi**2 / (4*self.height**2) )
                term_3 = math.sin( (2*j - 1)*math.pi*position / (2*self.height) )
                summationResult += term_1 * term_2 * term_3
            barP0 = self.alpha * self.Q / ( self.M + self.alpha * self.alpha * self.Q ) * \
                ( self.tao_0 + 0.5 * self.rho * self.g * self.height ) - self.rho_f * self.g * \
                ( 0.5 * self.height );
            displacementValue = 8.0*self.alpha*self.height*barP0*summationResult / (self.M*math.pi**2)
            displacementValue -= ( self.g / self.M ) * ( self.rho - self.alpha * self.rho_f ) * ( self.height * yPosition - 0.5 * yPosition**2.0 )
            displacementValue -= (self.tao_0*yPosition/self.M)
            return displacementValue


def pEquilibrium(stress, height, c_f, c_s, phi, ni, G, alpha):
    Lambda = 2*G*ni/(1-2*ni)
    Psi = phi*c_f + ( alpha - phi )*c_s
    Beta = Lambda*(1-ni)/ni + alpha*alpha/Psi
    p = -alpha*stress/( Psi*Beta )
    return p

def uEquilibrium(stress, height, c_f, c_s, phi, ni, G, alpha):
    Psi = phi*c_f + ( alpha - phi )*c_s
    lame = 2*G*ni/(1-2*ni)
    Beta = lame*(1-ni)/ni + alpha*alpha/Psi
    return stress*height/Beta
