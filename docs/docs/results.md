# Results

Understanding the energy distribution and number density of neutrons within the tank is important for the following two reasons:

1. <b> Activation Efficiency </b>: The effectiveness of the activation process depends on the energy of the neutrons. Relevant results are presented in sections [Flux](#1-flux) and [Energy Spectrum](#2-energy-spectrum).

2. <b> Safety Assurance </b>: Neutron energies at the tank boundary need to be low enough to ensure safe release into the environment. Findings related to this aspect are discussed in section [Escaping Neutrons](#3-escaping-neutrons).

Simulations involving 50000 neutrons were conducted to produce the results. Note that the neutron source emits approximately $8 \times 10^5$ neutrons per second. If applicable, the obtained result can be adjusted accordingly based on this emmision rate by transformation, or alternativaly, simulations utilizing the emmission rate as the number of neutrons could be conducted.

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

The flux as a function of distance from the origin was calculated using the methodology outline in section [Flux](./theory.md#61-flux) and is shown in Figure 14.  


<center>
<figure markdown="span">
  ![Flux](./figures/result_flux.png){ width="550" }
  <figcaption><i>Figure 14: Flux [number of neutrons per meter square] as a function of distance from the origin [meters]. </i></figcaption>
</figure>
</center>

## 2 Energy Spectrum

The energy spectrum is shown as a probability density graph for varying distances from the origin in Figure 15. 

<center>
<figure markdown="span">
  ![Flux](./figures/result_energy_spectrum.png){ width="550" }
  <figcaption><i>Figure 15: Probability density graph of energies at varying distances r from the origin. </i></figcaption>
</figure>
</center>

For neutron activation, we are not really interested in the thermal neutrons (left peak). In Figure 16 a movie of the energy spectrum is shown for varying distance $r$ from the origin, where a cut was made to highlight the thermal and high energy peak. 

<center>
<figure markdown="span">
  ![Flux](./figures/neutron_energy_spectrum.gif){ width="1000" }
  <figcaption><i>Figure 16: Movie of the probability density graph of energies for varying distances r from the origin. A cut was made to highlight the thermal and high energy peak. </i></figcaption>
</figure>
</center>

The energy spectrum and the [result for the flux](#1-flux) can be combined to find the number of neutrons that have a specific energy.

## 3 Escaping Neutrons

