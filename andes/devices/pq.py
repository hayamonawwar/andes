import logging
from collections import OrderedDict
from andes.core.model import Model, ModelData  # NOQA
from andes.core.param import IdxParam, NumParam # NOQA
from andes.core.var import Algeb, ExtAlgeb  # NOQA
from andes.core.discrete import Comparer  # NOQA
logger = logging.getLogger(__name__)


class PQData(ModelData):
    def __init__(self):
        super().__init__()
        self.bus = IdxParam(model='Bus', info="linked bus idx", mandatory=True)
        self.owner = IdxParam(model='Owner', info="owner idx")

        self.p0 = NumParam(default=0, info='active power load', power=True, tex_name=r'p_0')
        self.q0 = NumParam(default=0, info='reactive power load', power=True, tex_name=r'q_0')
        self.vmax = NumParam(default=1.1, info='max voltage before switching to impedance', tex_name=r'v_{max}')
        self.vmin = NumParam(default=0.9, info='min voltage before switching to impedance', tex_name=r'v_{min}')


class PQ(PQData, Model):
    def __init__(self, system=None, config=None):
        PQData.__init__(self)
        Model.__init__(self, system, config)
        self.group = 'StaticLoad'
        self.flags.update({'pflow': True})
        self.config.add(OrderedDict((('pq2z', 1), )))

        self.tex_names.update({'vcmp_zl': 'z_{vl}',
                               'vcmp_zi': 'z_{vi}',
                               'vcmp_zu': 'z_{vu}'})

        self.a = ExtAlgeb(model='Bus', src='a', indexer=self.bus, tex_name=r'\theta')
        self.v = ExtAlgeb(model='Bus', src='v', indexer=self.bus, tex_name=r'V')

        self.vcmp = Comparer(var=self.v, lower=self.vmin, upper=self.vmax,
                             enable=self.config.pq2z)

        self.a.e_str = "u * (p0 * vcmp_zi + \
                             p0 * vcmp_zl * (v ** 2 / vmin ** 2) + \
                             p0 * vcmp_zu * (v ** 2 / vmax ** 2))"

        self.v.e_str = "u * (q0 * vcmp_zi + \
                             q0 * vcmp_zl * (v ** 2 / vmin ** 2) + \
                             q0 * vcmp_zu * (v ** 2 / vmax ** 2))"

        # Experimental Zone Below
        # self.v_ref = Algeb(info="Voltage reference for PI")
        # self.kp = Service()
        # self.ki = Service()

        # self.kp.e_str = "1"
        # self.ki.e_str = "1"
        # self.pi = PIController(self.v, self.v_ref, self.kp, self.ki,
        #                        info='PI controller for voltage')