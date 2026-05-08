# Introduction


## Thermal simulator (T)
```python
sim = Simulator_T(eq_heat, time_control, outputs)
```

## Mechanical simulator (M)
```python
sim = Simulator_M(eq_mom, time_control, outputs, caverns)
```

## Thermo-mechanical simulator (TM)
```python
sim = Simulator_TM(eq_mom, eq_heat, time_control, outputs)
```

## Fully coupled TM simulator
```python
sim = Simulator_Full(eq_mom, eq_heat, time_control, outputs, caverns)
```





![Representation of the coupling scheme.](../../images/coupling_models.svg){#fig-couplings width="400%"}