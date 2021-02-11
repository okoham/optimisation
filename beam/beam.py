import math
from functools import cached_property, lru_cache


MATERIALS = {
    "AL7010": dict(
        E=71000, 
        nu=0.33, 
        ftu=515, 
        fty=455, 
        fcy=440, 
        fsu=295, 
        rho=2.82e-6, 
        cmat=5,        # material cost €/kg
        cprod=5        # machining cost €/kg
    ),
    "AL2198": dict(
        E=76000, 
        nu=0.33, 
        ftu=495, 
        fty=430, 
        fcy=415, 
        fsu=270, 
        rho=2.69e-6, 
        cmat=10,        # material cost €/kg
        cprod=10        # machining cost €/kg
    ),
    "TI64": dict(
        E=110000, 
        nu=0.33, 
        ftu=900, 
        fty=800, 
        fcy=780, 
        fsu=520, 
        rho=4.5e-6, 
        cmat=60,        # material cost €/kg
        cprod=60        # machining cost €/kg
    ),
}

# TODO: add additional parameter dstab - distance of lateral stabilisers

class Cantilever(object):
    """
    Cantilever beam, I section, clamped at x=0, point load at x=L
    """

    def __init__(self, matname, L, h, tw, blf, tlf, buf, tuf, dstab=None):
        self.L = L 
        self.h = h
        self.tw = tw
        self.blf = blf 
        self.tlf = tlf
        self.buf = buf 
        self.tuf = tuf
        self.matname = matname
        self.material = MATERIALS[matname]

        # default: no stabilisers; use full length
        self.dstab = dstab or self.L


    @cached_property
    def a_lf(self):
        """Return lower flange area"""
        return self.blf * self.tlf

    @cached_property
    def a_uf(self):
        """Return upper flange area"""
        return self.buf * self.tuf

    @cached_property
    def a_w(self):
        """Return web area"""
        return (self.h - self.tlf - self.tuf) * self.tw

    @cached_property
    def area(self):
        """Return cross section area"""
        return self.a_lf + self.a_uf + self.a_w

    @cached_property
    def cg_w(self):
        """Return z coodinate of web cg"""
        return (self.tlf + self.h - self.tuf)/2

    @cached_property
    def cg_uf(self):
        """Return z coordinate of upper flange cg"""
        return self.h - self.tuf/2

    @cached_property
    def cg_lf(self):
        """Return z coordinate of web cg"""
        return self.tlf/2

    @cached_property
    def cg(self):
        """Return z coordinate of section cg"""
        # sum ai zi / sum ai
        return (self.cg_uf*self.a_uf + self.cg_w*self.a_w + self.cg_lf*self.a_lf) / self.area

    @cached_property
    def iyy(self):
        """Return bending inertia (rel. to cg)"""
        # upper flange
        i = self.buf * self.tuf**3 / 12
        i += self.a_uf * (self.cg_uf - self.cg)**2
        # lower flange
        i += self.blf * self.tlf**3 / 12
        i += self.a_lf * (self.cg_lf - self.cg)**2
        # web 
        i += self.tw * (self.h - self.tuf - self.tlf)**3 / 12
        i += self.a_w * (self.cg_w - self.cg)**2
        return i

    @cached_property
    def it(self):
        """Return torsional inertia. Use thinwalled assumptions"""
        return (1/3) * (self.tw**3 * self.h + self.tlf**3 * self.blf + self.tuf**3 * self.buf)

    def mass(self):
        """Return beam mass"""
        return self.material['rho'] * self.L * self.area

    def cost(self):
        """Return beam cost"""
        # box we need to buy
        PAD = 3 # mm, pad each side
        mbox = self.material['rho'] * (self.L + 2*PAD) * (max(self.blf, self.buf) + 2*PAD) * (self.h + 2*PAD)
        # cost
        return mbox * self.material['cmat'] + (mbox - self.mass()) * self.material['cprod']

    #@lru_cache
    def Qz(self, F, x):
        """Return shear force at x"""
        if (x < 0) or (x > self.L):
            return math.nan        
        return F

    #@lru_cache
    def My(self, F, x):
        """return bending moment at x"""
        if (x < 0) or (x > self.L):
            return math.nan        
        return -F*(self.L - x)

    #@lru_cache
    def stress(self, F, x, z):
        """Return normal stress $\sigma_x$ at coordinate z"""
        if (z < 0) or (z > self.h):
            return math.nan
        return self.My(F, x) * (z - self.cg) / self.iyy

    def w_max(self, F):
        """Return max. deflection (at x = L)"""
        return F * self.L**3 / (3 * self.material['E']  * self.iyy)

    def tau_web(self, F, x):
        """Return average shear stress in web."""
        return abs(self.Qz(F, x)/self.a_w)
    
    def rf_t_uf(self, F):
        """Return reserve factor, upper flange tension"""
        smax = max([self.stress(F, 0, zi) for zi in (self.h, self.h - self.tuf)])
        if smax > 0: 
            return self.material['ftu'] / smax
        else:
            return math.inf

    def rf_c_uf(self, F):
        """Return reserve factor, upper flange compression"""
        smin = max([self.stress(F, 0, zi) for zi in (self.h, self.h - self.tuf)])
        if smin < 0: 
            return self.material['fcy'] / abs(smin)
        else:
            return math.inf

    def rf_lb_uf(self, F):
        """Return reserve factor, upper flange local buckling"""
        smin = min([self.stress(F, 0, zi) for zi in (self.h - self.tuf, self.h)])
        if smin >= 0:
            return math.inf
        else:
            b = (self.buf - self.tw)/2
            sallow = sigma_fb(b, self.tuf, self.material['E'])
            return abs(sallow/smin)        

    def rf_t_lf(self, F):
        """Return reserve factor, lower flange tension"""
        smax = max([self.stress(F, 0, zi) for zi in (0, self.tlf)])
        if smax > 0: 
            return self.material['ftu'] / smax
        else:
            return math.inf

    def rf_c_lf(self, F):
        """Return reserve factor, lower flange compression"""
        smin = max([self.stress(F, 0, zi) for zi in (0, self.tlf)])
        if smin < 0: 
            return self.material['fcy'] / abs(smin) 
        else:
            return math.inf

    def rf_lb_lf(self, F):
        """Return reserve factor, lower flange local buckling"""
        smin = min([self.stress(F, 0, zi) for zi in (0, self.tlf)])
        if smin >= 0:
            return math.inf
        else:
            b = (self.blf - self.tw)/2
            sallow = sigma_fb(b, self.tlf, self.material['E'])
            return abs(sallow/smin)

    def rf_s_web(self, F):
        """Return reserve factor, web material strength (von mises)"""
        tau = 1.5*self.tau_web(F, 0)
        return self.material['fsu'] / tau


    def rf_wb(self, F):
        """Return reserve factor, web buckling (combined shear and bending)"""
        b = self.h - self.tuf - self.tlf
        # shear
        tau = self.tau_web(F, 0)
        tau_cr = tau_crit(b, self.tw, self.material['E'], self.material['nu'])
        rs = tau / tau_cr
        # bending
        supper = self.stress(F, 0, self.h - self.tuf)
        slower = self.stress(F, 0, self.tlf)
        sx = min(supper, slower)
        # for some very degenerate beams (thick flanges), the min bending stress
        # can be > 0. In that case, set rb to zero (no positive influence) 
        if sx >= 0:
            rb = 0
        else:
            scrit = sigmab_crit(b, self.tw, self.material['E'])
            rb = abs(sx/scrit)
        # combined: HSB 45113-01, Eq. 3-1
        return 1 / math.sqrt(rs**2 + rb**2)


    def rf_lat(self, F):
        """Return reserve factor, lateral stability (Klein, Leichtbau, fig. 18.13)"""
        
        E = self.material['E'] 
        G = E / (2 * (1 + self.material['nu']))
    
        fcrit = (4.2/self.dstab**2) * math.sqrt(E * G * self.iyy * self.it)
        return abs(fcrit/F)


    def analyse(self, loads):
        """
        Return a dictionary containing inputs and analysis results
        """

        d = {
            'wmax': 0,
            'rf_t_uf': math.inf,
            'rf_t_lf': math.inf,
            'rf_c_uf': math.inf,
            'rf_c_lf': math.inf,
            'rf_lb_uf': math.inf,
            'rf_lb_lf': math.inf,
            'rf_s_web': math.inf,
            'rf_wb': math.inf,
            'rf_lat': math.inf,
            'Fmax': 0
        }

        for F in loads:
            d['wmax'] = max(d['wmax'], abs(self.w_max(F)))
            d['rf_t_uf'] = min(d['rf_t_uf'], self.rf_t_uf(F))
            d['rf_t_lf'] = min(d['rf_t_lf'], self.rf_t_lf(F))
            d['rf_c_uf'] = min(d['rf_c_uf'], self.rf_c_uf(F))
            d['rf_c_lf'] = min(d['rf_c_lf'], self.rf_c_lf(F))
            d['rf_lb_uf'] = min(d['rf_lb_uf'], self.rf_lb_uf(F))
            d['rf_lb_lf'] = min(d['rf_lb_lf'], self.rf_lb_lf(F))
            d['rf_s_web'] = min(d['rf_s_web'], self.rf_s_web(F))
            d['rf_wb'] = min(d['rf_wb'], self.rf_wb(F))
            d['rf_lat'] = min(d['rf_lat'], self.rf_lat(F))
            d['Fmax'] = max(d['Fmax'], abs(F))

        d['mass'] = self.mass()
        d['cost'] = self.cost()
        d['area'] = self.area

        for attr in ['L', 'h', 'tw', 'blf', 'tlf', 'buf', 'tuf', 'matname']:
            d[attr] = getattr(self, attr)

        return d

    def analyse_single(self, F):
        """
        Return a dictionary results for single load case
        """

        d = {}
        d['wmax'] = self.w_max(F)
        d['rf_t_uf'] =self.rf_t_uf(F)
        d['rf_t_lf'] = self.rf_t_lf(F)
        d['rf_c_uf'] = self.rf_c_uf(F)
        d['rf_c_lf'] = self.rf_c_lf(F)
        d['rf_lb_uf'] = self.rf_lb_uf(F)
        d['rf_lb_lf'] = self.rf_lb_lf(F)
        d['rf_s_web'] = self.rf_s_web(F)
        d['rf_wb'] = self.rf_wb(F)
        d['rf_lat'] = self.rf_lat(F)

        d['mass'] = self.mass()
        d['cost'] = self.cost()

        return d


def sigma_von_mises(sx, tau):
    """Return von Mises equivalent stress for 2d plane stress state"""
    return math.sqrt(sx**2 + 3*tau**2)


def tau_crit(b, t, E, nu):
    """
    Return shear buckling stress for simply supported plate, AR = infinity, 
    or b/a = 0.
    HSB 45112-01, figure 1, curve 3 -> k = 5.3 (all edges simply supported)
    """
    ks = 5.3
    return ks * math.pi**2 * (t/b)**2 * E / (12 * (1 - nu**2))


def sigmab_crit(b, t, E):
    """
    Return bending buckling stress for simply supported plate, AR = infinity, 
    or b/a = 0. Loaded in pure bending, i.e. sigma_max = sigma_min
    HSB 45111-01, figure 1, curve 2 -> k* = 21.58 (all edges simply supported)
    """
    k = 21.58
    return k * E * (t/b)**2 


def sigma_fb(b, t, E):
    """flange buckling stress. Plate, 3 sides simply supported.
    HSB 54111-01, figure 2, case 14, $\alpha -> \inf$"""
    return 0.367 * E * (t/b)**2


if __name__ == '__main__':
    import itertools
    import pandas as pd
    import random

    hset = (60, 80, 100, 120, 140, 160, 180, 200)
    tset_w = (1, 2, 1, 4, 5, 6)
    tset_f = (0, 3, 6, 9, 12)
    bset = (20, 30, 40, 50, 60, 70, 80)
    mset = ('AL7010', 'AL2198', 'TI64')

    L = 2000
    loads = [6000, -10000]

    # approach 1: all combinations
    if 1==0:
        all_res = []
        for (h, tw, blf, tlf, buf, tuf, matname) in itertools.product(hset, tset_w, 
                                                                    bset, tset_f, 
                                                                    bset, tset_f, mset):
            beam = Cantilever(matname, L, h, tw, blf, tlf, buf, tuf) 
            res = beam.analyse(loads)                                                                  
            all_res.append(res)
        
        df = pd.DataFrame(all_res)
        df.to_csv("beam_res_1.csv", index=False)
    
    # approach 2: random rampling
    N = 100000
    all_res = []
    for i in range(N):
        h = random.uniform(10, 300)
        tw = random.uniform(1.5, 6)
        blf = random.uniform(6, 100)
        tlf = random.uniform(0, 12)
        buf = random.uniform(6, 100)
        tuf = random.uniform(0, 12)
        matname = random.choice(mset)

        beam = Cantilever(matname, L, h, tw, blf, tlf, buf, tuf) 
        res = beam.analyse(loads)                                                                  
        all_res.append(res)

    df = pd.DataFrame(all_res)
    df.to_csv("beam_res_2.csv", index=False)
    