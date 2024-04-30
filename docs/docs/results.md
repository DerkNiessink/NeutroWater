# Results

Understanding the energy distribution and number density of neutrons within the tank is important for the following two reasons:

1. <b> Activation Efficiency </b>: The effectiveness of the activation process depends on the energy of the neutrons. Relevant results are presented in sections [Flux](#1-flux) and [Energy Spectrum](#2-energy-spectrum).

2. <b> Safety Assurance </b>: Neutron energies at the tank boundary need to be low enough to ensure safe release into the environment. Findings related to this aspect are discussed in section [Escaping Neutrons](#3-escaping-neutrons).

Simulations involving 50000 neutrons were conducted to produce the results. Considering that the neutron source emits approximately $8 \times 10^5$ neutrons per second, the outcomes were adjusted accordingly for ease of comparison. All simulation parameters that were used are shown below:

``` 
from neutrowater import diffusing_neutrons as dn

params = dn.Parameters(
            nNeutrons=50000, 
            radius_tank=0.225, 
            height_tank=0.85, 
            position_tank=(0, 0, -0.175)
            )
``` 

More information on how to reproduce these results is discussed in the [User Guide](user_guide.md#user-guide).

## 1 Flux

## 2 Energy Spectrum

## 3 Escaping Neutrons