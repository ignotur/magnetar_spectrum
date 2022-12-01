Magpies code 
====

MAGnetar Polarisation, Intensity and Spectroscopy code


## Compilation

To compile the code run:

```
make
```

## Procedure

The magnetosphere has to be solved first. In order to solve the equation for the twisted magnetosphere run:

```
python3 solve_twisted_magnetosphere.py 20 0.97 
```

where the first parameter is the degree of trigonometric polynomial used to expand the solution and the second parameter is p.

## Usage

The code should produce a photon list with phases. 
Current (simplified) version of the code can be ran as:
```
./magpies.out 
```
This simple run should create file `photon_list.txt` which contains the following: (1) theta, (2) phi and (3) final energy of photon only for photons which leave the system and reaches the observer

## References

Fernandez & Thompson (2007), ApJ, 660, 615

Nobili, Turolla & Zane (2008), MNRAS, 386, 1527

Pavan, Turolla, Zane & Nobili (2009), MNRAS, 395, 753

Thompson, Lyutikov & Kulkarni (2002), ApJ, 574, 332
