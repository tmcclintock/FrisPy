"""
This file creates a 'model object' which is initalized with 
a Frisbee's parameter values. This file calculates the
amplitudes of the lift force, drag force, and
torques about each principle axis.

For more information see, e.g. Hummel 2003.
"""
import numpy as np
"""
Constants:
PI
alpha_0: angle of attack of minimum drag. See Hummel 2003
"""
PI       = np.pi
alpha_0 = 4.0*PI/180. 

class Model(object):
	def __init__(self,*args):
                """
                Constructor

                PL0: lift parameter at 0 angle of attack (alpha)
                PLa: linear lift parameter that multiplies angle of attack
                PD0: drag parameter at alpha_0
                PDa: quadratic drag parameter, multiplies square angle of attack
                PTxwx: x-axis torque parameter linear in x-axis angular velocity
                PTxwz: x-axis torque parameter linear in z-axis angular velocity
                PTy0: y-axis torque parameter (pitch) at alpha = 0
                PTya: y-axis torque parameter linear in alpha
                PTywy: y-axis torque parameter linear in y-axis angular velocity
                PTzwz: z-axis torque parameter linear in z-axis angular velocity
                """
                if isinstance(args[0],np.ndarray): args = args[0]
                PL0,PLa,PD0,PDa,PTxwx,PTxwz,PTy0,PTya,PTywy,PTzwz=args
		self.PL0=PL0 
		self.PLa=PLa 
		self.PD0=PD0 
		self.PDa=PDa 
		self.PTxwx=PTxwx 
		self.PTxwz=PTxwz 
		self.PTy0=PTy0 
		self.PTya=PTya 
		self.PTywy=PTywy 
		self.PTzwz=PTzwz 

	def get_model_as_array(self):
		"""Return the model parameters just as an array of numbers.

		Args:
		    None
		    
		Returns:
		    params (array_like): Array of all model parameters.

		"""
		return np.array([self.PL0, self.PLa, self.PD0, self.PDa,
				 self.PTxwx, self.PTxwz, self.PTy0,
				 self.PTya, self.PTywy, self.PTzwz])

        def __str__(self):
                outstr = "Model:\n"+\
                         "PL0=%.2e\tPLa=%.2e\n"%(self.PL0,self.PLa)+\
                         "PD0=%.2e\tPDa=%.2e\n"%(self.PD0,self.PDa)+\
                         "PTxwx=%2.e\tPTxwz=%.2e\n"%(self.PTxwx,self.PTxwz)+\
                         "PTy0=%.2e\tPTwa=%.2e\tPTywy=%.2e\n"%(self.PTy0,self.PTya,self.PTywy)+\
                         "PTzwz=%.2e"%(self.PTzwz)
                return outstr

        """
        Coefficient functions.
        """
	def C_lift(self,alpha):
		return self.PL0 + self.PLa*alpha

	def C_drag(self,alpha):
		return self.PD0 + self.PDa*(alpha-alpha_0)*(alpha-alpha_0)

	def C_x(self,wx,wz):
		return self.PTxwz*wz + self.PTxwx*wx

	def C_y(self,alpha,wy):
		return self.PTy0 + self.PTywy*wy + self.PTya*alpha

	def C_z(self,wz):
		return self.PTzwz*wz

#An example of initializing and printing a model
if __name__ == "__main__":
        test_model = Model(0.33,1.9,0.18,0.69,0.43,-1.4e-2,-8.2e-2,-1.2e-2,-1.7e-3,-3.4e-5)
        print test_model
