# 2 & 3 body kinematic calculator for reactions in inverse kinematics
This project was developed to study (d,n) reactions in inverse kinematics for the RESONEUT detector array located at Florida State University's Fox Laboratory. 

## Project layout
```
├── build
│   ├── CMakeFiles
│   └── sim (executable)
├── CMakeLists.txt
├── data
│   └── test.root
├── files
│   └── nuclear_masses.json
├── include
│   ├── global.hpp
│   ├── kinematics.hpp
│   └── simulation.hpp
├── input
│   ├── dwba_l0.dat
│   └── test_input.json
├── README.md
├── source
│   ├── kinematics.cpp
│   └── simulation.cpp
├── tools
│   └── make_mass_table.py
└── vendor
    └── jsoncpp
        ├── json
        │   ├── json-forwards.h
        │   └── json.h
        └── jsoncpp.cpp
        └── LICENSE
```

## How to build & run this program...
```
mkdir build && cd build && cmake ../ && cd ..
```

```
cmake --build ./build/ --config Release -j(# of cores)
```

```
./build/sim --help
./build/sim input/test_input.json
```

## Future plans...
* Add s user input flag to run the efficiency 
* Add the silicon telescope geometry in input file
* Add a user input flag to check for particles incident on the detectors



