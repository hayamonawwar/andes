# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 12:52:21 2023

@author: haya.monawwar
"""

from andes.core import (Algeb, ConstService, ExtParam,
                        IdxParam, Lag, LeadLagLimit, Model, ModelData, NumParam, Washout)
import math


class TCSCBaseData(ModelData):
    """
    Base data for thyristor controlled series compensators.
    """

    def __init__(self):
        super().__init__()
        self.tcsc = IdxParam(model='TCSC',
                             info='Thyristor controlled series compensator idx',
                             mandatory=True,
                             unique=True,
                             )
        self.Sn = NumParam(default=100.0,
                           info="Power rating",
                           tex_name='S_n',
                           unit='MVA',
                           )
        self.Vn1 = NumParam(default=110.0,
                            info="AC voltage rating",
                            non_zero=True,
                            tex_name=r'V_{n1}',
                            unit='kV',
                            )
        self.Vn2 = NumParam(default=110.0,
                            info="rated voltage of bus2",
                            non_zero=True,
                            tex_name=r'V_{n2}',
                            unit='kV',
                            )

        self.P = NumParam(info='TCSC power rating. Equal to `Sn` if not provided.',
                          tex_name='P',
                          unit='MVA',
                          default=0.0,
                          )
        self.wref0 = NumParam(info='Base speed reference',
                              tex_name=r'\omega_{ref0}',
                              default=1.0,
                              unit='p.u.',
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
        self.wref = Algeb(info='Speed reference variable',
                          tex_name=r'\omega_{ref}',
                          v_str='wref0',
                          e_str='wref0 - wref',
                          )

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
                               tex_name=r'{\alpha}min',
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
                               tex_name='K{x}',
                               v_str='math.sqrt(Xc/Xl)'
                               )

        # FIXME: R should not be a constant because it depends on \alpha
        self.R = ConstService(info='Reactance of TCSC',
                              tex_name=r'X_TCSC(\alpha)',
                              v_str='[Xc * math.cos(Kx*(math.pi - alpha) \
                           * [(math.pi() * pow (Kx,4) - math.pi - 2*pow(Kx,4) \
                               *alpha + 2*pow(Kx,2)*alpha - math.sin(2*alpha)*pow(Kx,4) \
                            + math.sin(2*alpha)*pow(Kx,2) - \
                            4*pow(Kx,3)*pow(math.cos(alpha),2)*math.sin(math.pi() - alpha)\
                            - 4*pow(Kx,2)*math.cos(alpha)*math.sin(alpha) ))] ] / \
                              [math.pi * (pow(Kx,4) - 2*pow(Kx,2) + 1) \
                             * math.cos (Kx * (math.pi - alpha))]'
                              )

        self.pref = Algeb(info='Reference power input',
                          tex_name='P_{ref}',
                          v_str='pref',  # How to initialize this?
                          e_str='pref * K - pref',  # FIXME: equation error
                          )

        self.pbus = Algeb(info='Power at bus h',  # How to initialize this?
                          tex_name='P_{h}',
                          unit='p.u.',
                          )

        self.Washout = Washout(u='pref - pbus',
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
        self.output = Algeb(info='Final Output', tex_name=r'b(\alpha)',
                            e_str='1/R - output')


class TCSC(TCSC1Data, TCSC1Model):
    def __init__(self, system, config):
        TCSC1Data.__init__(self)
        TCSC1Model.__init__(self, system, config)
