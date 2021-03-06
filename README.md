# A Tighter Relaxation for the Relative Pose Problem Between Cameras

Convex relaxtions for the relative pose problem 
between two calibrated cameras. 
Refer to our paper [HERE](https://trebuchet.public.springernature.app/get_content/4be49395-75d2-4362-91d5-0fa53fce3660) 
for more information.



**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Jesus Briales](https://mapir.isa.uma.es/mapirwebsite/?p=2139), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GNUv3](https://github.com/mergarsal/TighterRPp/blob/main/LICENSE)


If you use this code for your research, please cite:

```
@ARTICLE{,
    author = {Garcia-Salguero, Mercedes and Briales, Jesus and Gonzalez-Jimenez, Javier},
     month = {April},
     title = {A Tighter Relaxation for the Relative Pose Problem Between Cameras},
   journal = {Journal of Mathematical Imaging and Vision},
    volume = {},
    number = {},
      year = {2022},
       url = {https://link.springer.com/article/10.1007/s10851-022-01085-z},
       doi = {10.1007/s10851-022-01085-z},
     pages = {}
}
```



# Dependencies

The solvers require *SDPA*. 

*SDPA*: 
    web: http://sdpa.sourceforge.net
    Download: http://sdpa.sourceforge.net/download.html#sdpa
    

## Build
```
git clone https://github.com/mergarsal/TighterRPp.git
cd TighterRPp

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. Run: 
```
        ./bin/example_base
```
 


# Launch experiments (tests)
cpulimit -l 45 ./../../build/bin/generic_test

*Note*: the experiments output several files. 
We suggest to run them on a different folder.
