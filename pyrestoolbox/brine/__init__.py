from .brine import *
from . import pr_vphi
from .pr_vphi import V2_inf as V_phi_pr
from . import vphi_route
from .vphi_route import V_phi, route_used

# Gas-free brine viscosity: IAPWS-2008 water x ion-additive Jones-Dole salt
# ratio x Kestin's measured pressure factor. Multi-salt over the pitzer.dat
# ions; a salinity given with no species named means NaCl, matching the
# density side. `route='mao_duan'` reproduces the pre-3.7.3 chain.
from . import viscosity_route
from .viscosity_route import brine_viscosity, salt_ratio
from . import iapws_viscosity
from .iapws_viscosity import mu_iapws2008
from . import appelo_volumes
from . import jones_dole_viscosity
from . import kestin_nacl_viscosity
