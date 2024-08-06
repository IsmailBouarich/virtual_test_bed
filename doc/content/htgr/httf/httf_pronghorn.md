# HTTF pronghorn

*Contact: Ismail Bouarich, Lise Charlot (Lise.Charlot.at.inl.gov), Mauricio Tano (Mauricio.Tano.at.inl.gov)*

## 3D porous media model description 

This is a 3D porous media model of Oregon State University's High Temperature Test Facility (HTTF), a one-quarter scale model of the General Atomics’ Modular High Temperature Gas-cooled Reactor (MHTGR). The HTTF uses electrical heaters as power source and helium as coolant and operates at prototypical temperature but at reduced pressure. The following figure shows two main components of the HTTF : the Reactor Pressure Vessel (RPV) and Reactor Cavity Cooling System (RCCS). The RCCS is not modeled in this work. The primary flow components in the HTTF include the top plenum, top reflector, core, and bottom reflector. In the radial direction, the model consists of layers starting from the core (including upper and lower reflectors), followed by the permanent reflector spanning from the top of the upper reflector to the bottom of the lower reflector. Axially, the model includes the helium gap and core barrel spanning across all these flow components.

!media /httf/httf_pronghorn/httf_cad_model.png
    style=width:50%
    caption=HTTF CAD model (source : [Oregon State University](https://engineering.oregonstate.edu/NSE/research-innovation/facilities/high-temperature-test-facility-httf))
    

The porous media model is a computational approach used to simulate the behaviour of the flow within the core. This model simplifies the representation of the complex core elements into porous media regions. Each region is characterized by its porosity, which is the fraction of void space, permeability, a measure of how easily fluid can flow through the medium, and hydraulic diameter. This allows for a considerable reduction in computational cost. The porous medium approach has been used before to model pebble bed reactors, and it can also be used to model prismatic high temperature reactors.

!media /httf/httf_pronghorn/homo.png
    style=width:40%
    caption=Regional homogenization of the HTTF core

!table id=ring_bc caption=Characteristics of each homogenized region
|  Region | Porosity | Hydraulic diameter (m)  |
| :- | :- | :- |
| 0 | 0 | 0 |
| 1 | 0.04047856 | 0.01905000 | 
| 2 | 0.02023928 | 0.00952500 | 
| 3 | 0.03822975 | 0.01079500 | 
| 4 | 0.16866067 | 0.01587500 | 
| 5 | 0.14842139 | 0.01496786 | 
| 6 | 0.07421070 | 0.01164167 | 
| 7 | 0.02023928 | 0.00952500 | 
| 8 | 0 | 0 |
| 9 | 0.02811011 | 0.01587500 | 

To solve the problem, the MOOSE Navier Stokes module is used. The top plenum and helium gap are assumed to operate under free flow conditions. Solid heat conduction is simulated for all components except the top plenum and the helium gap. The inlet plenum acts as an inlet flow distributor with a porosity of 1 (same as helium gap). 

The simulations of all models are done under steady state conditions which can be found in the table below.

!table id=ring_bc caption=Steady state operating conditions.
|  Parameter | Unit | Value  |
| :- | :- | :- |
| Heating power | $MW$ | 2.2 |
| Helium mass flow rate    | $kg/s$ | 1.0 |
| Helium inlet temperature    | $K$ | 500 |
| Helium pressure    | $MPa$ | 0.7 |


## Subscale model description    

The subscale model represents a detailed assembly for each hexagonal region in the core. Assemblies 2, 3, 4, 5, 6 and 7 were all modeled and included in this work. The coupling between the three-dimensional porous media model and the sub-assemblies allows to predict the heat distribution at a modest computational cost, compared to a fully meshed and detailed prototype. The assemblies solve the heat conduction, gap heat transfer and provide the wall temperature for the coolant channels. 

!media /httf/httf_pronghorn/ass2.png
    style=width:22%
    caption=Assembly 2

!media /httf/httf_pronghorn/ass3.png
    style=width:25%
    caption=Assembly 3

!media /httf/httf_pronghorn/ass4.png
    style=width:25%
    caption=Assembly 4    

!media /httf/httf_pronghorn/ass5.png
    style=width:25%
    caption=Assembly 5

!media /httf/httf_pronghorn/ass6.png
    style=width:25%
    caption=Assembly 6

!media /httf/httf_pronghorn/ass7.png
    style=width:22%
    caption=Assembly 7        

## HTTF Multi App model description

The coupling between the porous medium model and the subscale model was done thanks to the MultiApp system in MOOSE. This system allows the porous media main App to transfer solid temperature on the boundaries of the hexagonal regions to the different subscale models for each region, thus exploring in more detail the temperature distribution in each assembly. The sub-assemblies are themselves coupled to a 1-D channel flow sub-App that solves the helium flow in the coolant channels and transfers the values back to the assemblies. 

!media /httf/httf_pronghorn/Coupling_strategy.png
    style=width:30%
    caption=Coupling strategy

## Input File Description

### Initialization 

The first lines in the input file are general parameters governing the system. The inlet temperature, outlet pressure and power are imposed. The material properties are also given. The porosity and hydraulic diameter are calculated. The air to volume fraction is a number specific to each region which converts the wall heat transfer coefficient from a surface parameter to a volume parameter as required by the model. Follows the block categorization where flow blocks, solid blocks and power blocks are defined. 

### Navier Stokes

This block solves the Navier Stokes equation in the flow blocks using the finite volume method. The problem is considered weakly compressible and it is naturally treated as a porous medium. Gravity is also included but could be excluded as it doesn't considerably affect the soluition. The initial and boundary conditions are also specified. Darcy and Forchheimer are the considered friction types.  

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Modules language=cpp

### Variables

This block defines the global variables governing the problem which are the fluid temperature, and solid temperature.

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Variables language=cpp

### FVKernels

This section provides finite volume kernels used to solve for 3D heat conduction. A comment explains what each kernels adds in the block. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=FVKernels language=cpp

### FVBCs

This part provides finite volume boundary conditions which govern 3D heat conduction in the system. A comment explains each boundary condition. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=FVBCs language=cpp

### Functions

This includes the relevant spacial and time dependent functions for the system such as the inlet velocity, outlet pressure and power. Some function are included in the common.i file as they are common the MainApp and SubApps. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Functions language=cpp

### FunctorMaterials

Material properties are input in this block. A comment explains the feature implemented by each block. The first one for instance adds material properties (specific heat, thermal conductivity...l) to the system from fluid properties. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=FunctorMaterials language=cpp

### AuxVariables

This block creates relevant auxiliary variables to the problem, such as power density, porosity and interstitial velocity. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=AuxVariables language=cpp

### AuxKernels

The evolution of the auxiliary variables is governed by this block. A comment explains each auxKernel

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=AuxKernels language=cpp

### Postprocessors

This part outputs relevant data from the simulation to be further processed in Excel or Python.

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Postprocessors language=cpp

### Preconditioning

This block describes the preconditioner used by the solver. 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Preconditioning language=cpp

### Executioner 

This block describles the Newthon solver calculation process, where start time, end time and adaptative time stepper are specified. A nonlinear absolute tolerance of 1e-1 and relative tolerance of 1e-4 are imposed. Automatic scaling is 

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Executioner language=cpp

### MultiApps

This block is where subassemblies are connected to the 3D porous media model via the MOOSE MultiApp system. For assemblies 2, 3, 5, 6 and 7, a rotation is needed to insert each assembly accurately in the core. Notably, assemblies 2 and 7 are the same within 180 degrees rotation, which is why subscale 2 positions are included in subscale 7 SubApp. As the hexagon has 6 different faces, each subscale SubApp represents one assembly in one direction. For instance, subscale_3_1 represents the assemblies 3 which are all face the same direction, subscale_3_2 represents the assemblies 3 which are 60 degres rotated clockwise to those in subscale_3_1 and so on. The positions of each subscale assembly is either unique and directly implemented or is specified via a file which contains all the positions of the included assemblies.  

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=MultiApps language=cpp

### Transfers

In this block, transfers between the MainApp and the SubApps are defined. The outlet pressure and the mass flow rate are transfered from the MainApp to the SubApps, in order to have a coherent coupling. The boundary temperatures of all hexagons are transfered to the boundaries of the assemblies concerned. The temperatures on theses boundaries represent boundary conditions for the subscale model.    

!listing htgr/httf/httf_pronghorn/HTTF_steady_state.i block=Transfers language=cpp

## Subscale Input File Description 

Here, the input files of subscale models are described. As the files have a very similar baseline, only the assembly 4 file will be detailled.  

### Mesh

This is the main block for building the subscale model. Firstly, fuel rod pincells, a large coolant pincell and a dummy pincell are created using the PolygonConcentricCircleMeshGenerator. Dummy pincells were used to complete the shape needed to form a hexagonal assembly using the PatternedHexMeshGenerator. The HexagonMeshTrimmer then trims the peripheral boundary of the assembly, preparing it for a further modification. A rotation is needed for the model to match the HTTF core assembly 4, followed by an extrusion in the z direction. Then, dummy pincells are deleted and plane cuts are performed to create the desired assembly. The gaps between the heater rods and the ceramic block are built by pairing 2 subdomains for each gap (inner and outer parts of gap) and deleting the block between them. External boundaries of the assembly are defined in assembly_boundary block using SideSetsFromNormalsGenerator. 

!listing htgr/httf/httf_pronghorn/assembly_4.i block=Mesh language=cpp

### Heat conduction 

The temperature variable T is created and initialized. Additional auxiliary variables, such as fluid temperature, power density, boundary temperature and wall temperature are generated to fully describe the heat conduction in the assembly. Kernels and AuxKernels explicit the equations that govern these variables.  

### Boundary conditions 

There are 2 boundary conditions in the system. The first one couples far field T_fluid received by the coolant channels to the temperature of the solid in the assembly with a heater transfer coefficient htc. The second one couples the boundary temperature received from the porous media model (MainApp) to the external boundaries of the assembly. 

!listing htgr/httf/httf_pronghorn/assembly_4.i block=BCs language=cpp

### UserObjects

The NearestPointLayeredSideAverage userobject is used to sample in 50 layers the temperature in the wall of the coolant channels in the assembly MainApp in the z direction, which will then be transfered to the THM SubApp.   

!listing htgr/httf/httf_pronghorn/assembly_4.i block=UserObjects language=cpp

### Postprocessors

Multiple postprocessors are developped to provide more insight into all aspects of the heat conduction phenomenon in the assembly and verify the model.

!listing htgr/httf/httf_pronghorn/assembly_4.i block=Postprocessors language=cpp

### MultiApps

This block creates the 1D THM coolant channels SubApp. The positions of the coolant channels are specified in the large_coolant_position.txt. When coolant channels have different diameters, they are to be specified too. In this case, there are only large coolant channels. One processor is allocated for each coolant channel.  

!listing htgr/httf/httf_pronghorn/assembly_4.i block=MultiApps language=cpp

### Transfers 

In this part, transfers between the assembly and the coolant channels occur, thus coupling the two. The wall temperature is transferred from the MainApp to the 1D coolant channel, which gives back the coolant temperature and heat transfer coefficient after solving the flow in the SubApp. Mass flow rate and pressure are also transferred to the THM SubApp.  

!listing htgr/httf/httf_pronghorn/assembly_4.i block=Transfers language=cpp

## Coolant flow Input File Description 

The input file of the last part of the MultiApp pyramid coolant flow is explained here. Each coolant diameter has its own input file, however they have a very similar file design with changes in diameter and consequent parameters.

The code begins by specifying the physical parameters of the coolant channels, including the dimensions of different channel types, such as large, medium, and small channels. It calculates the flow areas and mass flow rates based on these dimensions.

The FluidProperties block specifies the characteristics of helium, which is considered here as an ideal gas with constant properties. The Closures section specifies the closure factors used which are the Churchill coefficients for wall friction factor and Dittus-Boelter coefficients for wall heat transfer coefficient. AuxVariables introduces the wall temperature, which is received from the assembly MainApp and the heat transfer coefficient, which is updated by the kernel htc_aux that calculates the average of the heat transfer property. The Components section outlines the different components of the simulation, including inlets, outlets, and the flow channel with specified direction, number of elements, area and hydraulic diameter, as well as the boundary condition of heat transfer from the wall temperature. The ControlLogic block allows for transfer and subsequent overwriting of mass flow rate and pressure from the Subscale MainApp. Postprocessors are computed to verify physical parameters of the model. 