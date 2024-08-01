![alt text](<./docs/source/_static/logo_2.png>){width=35%}

# Introduction
The SafeInCave simulator is developed for simulating gas storage operations in salt caverns. The finite element implementation uses the FEniCS package. Additionally, it includes different constitutive models for salt rock mechanics. 

## Getting started
The best way to get started is to visit the examples in the folder "apps" and run the "main.py" files. You can create your own examples by properly editing files "input_model.json" and "input_bc.json". In order to have a better idea of how the simulator works, check the function H2CaveSimulator at "libs/Simulators.py".

The simulator also includes a material point model, which is useful for calibration purposes. Check function SmpSimulator at "libs/Simulators.py".

## Current members : 
- [Hermínio Tasinafo Honório] (H.TasinafoHonorio@tudelft.nl),  Maintainer, 2023-present
- [Hadi Hajibeygi] (h.hajibeygi@tudelft.nl), Principle Investigator
