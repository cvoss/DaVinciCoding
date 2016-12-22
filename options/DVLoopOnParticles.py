##
# $Id: DVLoopOnParticles.py,v 1.1 2008-12-23 11:31:30 pkoppenb Exp $
#
# Run Loop on Particles algorithm
#
# @autor P. Koppenburg - J. Palacios
#
##
from Gaudi.Configuration import *

from Configurables import LoopOnParticles, PhysDesktop
loop = LoopOnParticles()
loop.HistoProduce = True
loop.addTool(PhysDesktop())
loop.PhysDesktop.InputLocations = [ "StdLoosePions" ]

from Configurables import DaVinci
DaVinci().EvtMax = 100                                                           # Number of events
DaVinci().DataType = "2008"                                                      # Default is "DC06"
DaVinci().Simulation   = True                                                    # It is MC
DaVinci().UserAlgorithms = [ loop ]                                              # Declare your aldo
DaVinci().HistogramFile = "DaVinciUser.hbook"
