# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 12:52:21 2023

@author: haya.monawwar
"""

from andes.core import (Algeb, ConstService, ExtParam,
                        IdxParam, Lag, LeadLagLimit, Model, ModelData, NumParam, Washout)

from andes.core.var import ExtAlgeb


class TCSCBaseData(ModelData):
    """
    Base data for thyristor controlled series compensators.
    """

    def __init__(self):
        super().__init__()
        self.bus1 = IdxParam(model='Bus', info="idx of from bus")
        self.bus2 = IdxParam(model='Bus', info="idx of to bus")

        self.Sn = NumParam(default=100.0,
                           info="Power rating",
                           tex_name='S_n',
                           unit='MVA',
                           )

        self.Vn1 = NumParam(default=110.0,
                            info="AC voltage rating",
                            non_zero=True,
                            tex_name=r'V_{h}',
                            unit='kV',
                            )

        self.Vn2 = NumParam(default=110.0,
                            info="AC voltage rating",
                            non_zero=True,
                            tex_name=r'V_{k}',
                            unit='kV',
                            )

        self.P = NumParam(info='TCSC power rating. Equal to `Sn` if not provided.',
                          tex_name='P',
                          unit='MVA',
                          default=0.0,
                          )

        self.alpha = NumParam(info='Firing angle',
                              tex_name='alpha',
                              default=0.0,
                              unit='p.u.',
                              )


class TCSCBase(Model):
    """
    Base TCSC model.

    """

    def __init__(self, system, config):
        Model.__init__(self, system, config)
        self.flags.update({'tds': True})

        self.Xc = NumParam(info='Capacitive reactance ',
                           tex_name='Xc',
                           unit='p.u.',
                           power=True,
                           )
        self.Xl = NumParam(info='Inductive reactance ',
                           tex_name='Xl',
                           unit='p.u.',
                           power=True,
                           )


class TCSC1Data(TCSCBaseData):
    def __init__(self):
        super().__init__()
        self.Kw = NumParam(info='Regulator gain ',
                           tex_name='K_w',
                           default=0.05,
                           unit='p.u.',
                           ipower=True,
                           )

        self.T1 = NumParam(info='Low-pass time constant',
                           default=0.1,
                           tex_name='T_1')
        self.T2 = NumParam(info='Lead time constant',
                           default=0.2,
                           tex_name='T_2')
        self.T3 = NumParam(info='Lag time constant',
                           default=10.0,
                           tex_name='T_3')
        self.Tw = NumParam(info='Washout time constant',
                           default=10.0,
                           tex_name='T_w')

        self.alphaL = NumParam(info='Maximum firing angle',
                               tex_name=r'{\alpha}_{max}',
                               default=0.0,
                               unit='p.u.',
                               )
        self.alphaU = NumParam(info='Minimum firing angle',
                               tex_name=r'{\alpha}_min',
                               default=0.0,
                               unit='p.u.',
                               )


class TCSC1Model(TCSCBase):
    """
    Implement TCSC1 model.
    """

    def __init__(self, system, config):
        TCSCBase.__init__(self, system, config)

        self.Kx = ConstService(info='Constant to calculate reactance of TCSC at a firing angle',
                               tex_name='K_{x}',
                               v_str='sqrt(Xc/Xl)'
                               )

        self.a1 = ExtAlgeb(model='Bus', src='a', indexer=self.bus1, tex_name='a_1',
                           info='phase angle of the from bus',
                           ename='Pij',
                           tex_ename='P_{ij}',
                           )
        self.a2 = ExtAlgeb(model='Bus', src='a', indexer=self.bus2, tex_name='a_2',
                           info='phase angle of the to bus',
                           ename='Pji',
                           tex_ename='P_{ji}',
                           )
        self.v1 = ExtAlgeb(model='Bus', src='v', indexer=self.bus1, tex_name='v_1',
                           info='voltage magnitude of the from bus',
                           ename='Qij',
                           tex_ename='Q_{ij}',
                           )
        self.v2 = ExtAlgeb(model='Bus', src='v', indexer=self.bus2, tex_name='v_2',
                           info='voltage magnitude of the to bus',
                           ename='Qji',
                           tex_ename='Q_{ji}',
                           )

        # FIXME: R should not be a constant because it depends on \alpha    // DONE.

        self.R = Algeb(info='Reactance of TCSC',
                       tex_name=r'X_TCSC(\alpha)',
                       v_str='[Xc * cos(Kx*(pi - alpha) '
                             '* [(pi * pow(Kx,4) -  pi - 2*pow(Kx,4) '
                             '*alpha + 2*pow(Kx,2)*alpha - sin(2*alpha)*pow(Kx,4)'
                             '+ sin(2*alpha)*pow(Kx,2) - '
                             '4*pow(Kx,3)*pow(cos(alpha),2)*sin(pi - alpha)'
                             '- 4*pow(Kx,2)*cos(alpha)*sin(alpha) ))] ] / '
                             ' [pi * (pow(Kx,4) - 2*pow(Kx,2) + 1) '
                             '* cos (Kx * (pi - alpha))]',
                       #    e_str='[Xc * cos(Kx*(pi - alpha) \
                       #        * [(pi * pow(Kx,4) - pi - 2*pow(Kx,4) \
                       #            *alpha + 2*pow(Kx,2)*alpha - sin(2*alpha)*pow(Kx,4) \
                       #         + sin(2*alpha)*pow(Kx,2) - \
                       #         4*pow(Kx,3)*pow(cos(alpha),2)*sin(pi - alpha)\
                       #         - 4*pow(Kx,2)*cos(alpha)*sin(alpha) ))] ] / \
                       #           [pi * (pow(Kx,4) - 2*pow(Kx,2) + 1) \
                       #          * cos (Kx * (pi - alpha))] - R',
                       )

        self.output = Algeb(info='Final Output', tex_name=r'b(\alpha)',
                            v_str='1/R',
                            e_str='1/R - output')

        # FIXME: P and Q at bus h and bus k need to be defined/initialized & make Pref const. // DONE.
        self.pref = ConstService(info='Reference power input',
                                 tex_name='P_{ref}',
                                 v_str='pref',  # How to initialize this?
                                 )

        self.Ph = Algeb(info='Active power at bus h',
                        tex_name='P_{h}',
                        unit='p.u.',
                        v_str='v1*v2*output*sin(a1 - a2)',
                        e_str='v1*v2*output*sin(a1 - a2) - Ph'
                        )

        self.Qh = Algeb(info='Reactive power at bus h',
                        tex_name='Q_{h}',
                        unit='p.u.',
                        v_str='pow(v1,2)*output - v1*v2*output*cos(a1 - a2)',
                        e_str='pow(v1,2)*output - v1*v2*output*cos(a1 - a2) - Qh'
                        )

        self.Pk = Algeb(info='Active power at bus k',
                        tex_name='P_{h}',
                        unit='p.u.',
                        v_str='-v1*v2*output*sin(a1 - a2)',
                        e_str='-v1*v2*output*sin(a1 - a2) - Pk'
                        )

        self.Qk = Algeb(info='Reactive power at bus k',
                        tex_name='Q_{h}',
                        unit='p.u.',
                        v_str='pow(v2,2)*output - v1*v2*output*cos(a1 - a2)',
                        e_str='pow(v2,2)*output - v1*v2*output*cos(a1 - a2) - Qk'
                        )

        self.Washout = Washout(u='pref - Ph',
                               K='Kw * Tw',
                               T=self.Tw,
                               tex_name='x1'
                               )

        self.Lag = Lag(u=self.Washout_y,
                       T=self.T1,
                       K=1, D=1,
                       tex_name='x2'
                       )

        self.LL = LeadLagLimit(u=self.Lag_y,
                               T1=self.T2,
                               T2=self.T3,
                               lower=self.alphaL,
                               upper=self.alphaU,
                               tex_name='x3'
                               )


class TCSC(TCSC1Data, TCSC1Model):
    def __init__(self, system, config):
        TCSC1Data.__init__(self)
        TCSC1Model.__init__(self, system, config)
