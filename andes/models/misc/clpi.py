
from andes.core import (ConstService, ExtParam, IdxParam, Model, ModelData, NumParam)

from andes.core.var import ExtAlgeb
from andes.core.block import PIDController




class CLPIData(ModelData):
    """
    Data for open-loop PI controller..
    """

    def __init__(self):
        super().__init__()
        self.gov = IdxParam(model='TurbineGov',
                            info='Turbine governor idx',
                            mandatory=True,
                            )
        self.kP = NumParam(info='PID proportional coeff.',
                           tex_name='k_P',
                           default=1,
                           )
        self.kI = NumParam(info='PID integrative coeff.',
                           tex_name='k_I',
                           default=1,
                           )
        self.kD = NumParam(info='PID derivative coeff.',
                           tex_name='k_D',
                           default=0,
                           )
        self.tD = NumParam(info='PID derivative time constant coeff.',
                           tex_name='t_D',
                           default=0,
                           )

class CLPIModel(Model):
    """
    Implementation for close-loop PI controller.
    """

    def __init__(self, system, config):
        Model.__init__(self, system, config)
        self.group = 'Experimental'
        self.flags.tds = True

        self.wd = ExtAlgeb(model='TurbineGov', src='pout', indexer=self.gov,
                           info='Generator speed deviation',
                           unit='p.u.',
                           tex_name=r'\omega_{dev}',
                           )
        self.pout = ExtAlgeb(model='TurbineGov', src='pout', indexer=self.gov,
                             tex_name='P_{out}',
                             info='Turbine governor output',
                             )
        self.pout0 = ConstService(v_str='pout',
                                  tex_name='P_{out0}',
                                  info='initial turbine governor output',
                                  )
        self.PID = PIDController(u=self.wd, kp=self.kP, ki=self.kI,
                                 kd=self.kD, Td=self.tD,
                                 tex_name='PID', info='PID', name='PID',
                                 ref=self.pout0,
                                 )
        self.pref = ExtAlgeb(indexer=self.gov,
                             tex_name='P_{ref}',
                             info='Turbine governor output',
                             model='TurbineGov',
                             src='Pref',
                             e_str='PID_y',
                             v_str='pref',
                             )

class CLPI(CLPIData, CLPIModel):
    r"""
    Closed-loop PI controller that takes Generator speed deviation as input.

    ```
        ┌────────────────────┐
        │      ki     skd    │
    u ->│kp + ─── + ───────  │ -> y
        │      s    1 + sTd  │
        └────────────────────┘
    ```
    """

    def __init__(self, system, config):
        CLPIData.__init__(self)
        CLPIModel.__init__(self, system, config)