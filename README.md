# PYTOP

# Summary 
## This python module calculates the topological invarint of 2 dimensional quantum hall systems in a computationally efficient way. Addionaly it also gives the edge states. Bulk boundary correspondence can be varified. 

### Installation
```bash
git clone https://github.com/kumarjnc/PYTOP.git
```

# Example
## Topological invariant calculation
```python
import qsh
import link
import field
import params
import numpy as np
```
```python
energy, energy_val,nth_m , nky_m=qsh.bandstr(params.tx,params.ty,params.Nk,params.Ntheta,params.size_bzone,params.alpha,params.lmda,params.gama,square=True)
```

```python
gap=np.min(abs(energy))
Nth, Nky, no_of_band_filled=energy.shape
no_of_band=int(no_of_band_filled/2)
evec=energy_val
```
```python
U1, U2=link.link(params.Nk, params.Ntheta, no_of_band, evec)
```
```python
mu=field.field(Nky, Nth, U1,U2)
if (gap > 0.01):
    muf=mu
else:
    muf=0.0
print('Z2 invariant %5.3f.' % muf)
```

