title:  Tutorial
author: Jabir Ali Ouassou
date:   2018-09-03

### Introduction
In this tutorial, we demonstrate how to model some common superconducting spintronics devices using the GNEUS software suite.
The first few sections will gradually introduce basic calculations, selfconsistent calculations, and nonequilibrium calculations.
After that, we will demonstrate some more advanced features of the configuration files, such as how to initialize a system using mathematical equations, and how to set system parameters using command-line arguments.
Finally, we will show how to use more specialized simulation programs to calculate phase diagrams and the superconducting critical temperature.
The instructions assume that you use Bash as your terminal shell, and have Gnuplot available for data visualization.

### Basic calculations
As the first example, we consider an S/N/S Josephson junction.
This consists of a normal-metal layer in-between two bulk superconductors, which we treat as reservoirs with fixed order parameters \(\Delta = \Delta_0 e^{\pm i\pi/4}\).
This produces a phase-difference \(\delta\varphi = \pi/2\) and maximizes the charge supercurrent \(J_{\mathrm{e}} = J_{\mathrm{e}0} \sin(\delta\varphi)\).
The normal-metal layer is assumed to have a length \(L=3\xi\), where \(\xi\) is the superconducting coherence length.
Finally, the interfaces are assumed to have tunneling conductances \(G_{\mathrm{T}} = 0.3G_{\mathrm{N}}\), where \(G_{\mathrm{N}}\) is the normal-state conductance of the normal-metal layer.

First, we have to write an INI-like configuration file that describes the physical system above.
The format used for these configuration files is simple but flexible: 

  * Sections like `[material]` defines a new material layer in a one-dimensional nanostructure.
    Possible materials include `superconductor`, `conductor`, `ferromagnet`, and `halfmetal`.
  * The physical properties and numerical model used for each material is configured using key-value pairs with the syntax `key: value`.
    Depending on the property, this can be either a boolean variable (T or F), an integer (e.g. 1), a real number (e.g. 1.0), or a real vector (e.g. [1.0,0,0]).
  * Material properties can also be specified using mathematical functions or command-line options.
    This functionality will be discussed in more detail later in the tutorial.
  * Spaces are ignored by the configuration parser, and empty lines are permitted as well.
    However, tabs are *not* permitted anywhere, and newlines are *not* permitted within a key-value pair.
  * Everything after `#` is regarded as a comment, and ignored by the configuration parser.

Here is a configuration file that describes the physical setup above:

    :::ini
    # S/N/S Josephson junction
    # with π/2 phase difference

    [superconductor]
      order:         0
      gap:           1.00
      phase:        +0.25
    
    [conductor]
      order:         1
      length:        3.00
      conductance_a: 0.30
      conductance_b: 0.30
    
    [superconductor]
      order:         0
      gap:           1.00
      phase:        -0.25

Sections like `[superconductor]` and `[conductor]` specify that we wish to simulate an S/N/S junction.
Below each section, the physical properties of that material layer is specified using key-value pairs.

Let us first consider their physical properties.
For the superconductors, the parameter `phase` sets the phases of the order parameters (in units of \(\pi\)), while the parameter `gap` sets the magnitude (in units of \(\Delta_0\)).
Since \(\Delta = \Delta_0\) is the default value for superconducting layers, specifying the `gap` is in this case optional.

For the normal metal, the parameters `conductance_a` and `conductance_b` define the tunneling conductance of the "left" and "right" interfaces, respectively.
Here, "left" refers to the material defined above, and "right" the material defined below.
These conductance values are normalized relative to the bulk conductance of the normal metal.
The length of the material is obviously defined using the `length` parameter (in units of \(\xi\)).

In addition to these physical parameters, a property `order` is used to control the numerical procedure.
The simulation programs work by solving the Usadel diffusion equation in one layer at a time, and this property can be used to control in what order this happens.
This is particularly useful in large multilayer structures, where some iteration orders may converge faster than others due to e.g. symmetry reasons.
The special value 0 is used to disable simulations entirely in a material, thus making it into a reservoir.
This is what we have done in the example above: since both superconductors have the order set to 0, they are treated as bulk reservoirs by the simulation program.
Note that since the order defaults to 1, it is optional to define it for the normal metal here.

After saving the above configuration file in a suitable folder as e.g.  `josephson.conf`, open a terminal in that directory.
The simulation can now be started using the `converge` program:

    :::bash
    converge josephson.conf &

The simulation program should now run in the background.
It will continuously write status information to the file `output.log`, while information about any errors that might occur is written to `error.log`.
In order to track the progress of the simulation, you can use the command `tail`:

    :::bash
    tail -n +0 -f *.log

If you compiled the program using IFort, and have a modern processor (e.g. an Intel Core i7 in my case), the simulation should finish in less than a minute.
During that time, it solves the Usadel equation as a function of position for 1000 energies in the range \((0,30\Delta_0)\).
Physical observables such as the density of states and supercurrents are subsequently calculated from the propagators, and the results saved to DAT files.

For this system, we may be especially interested in displaying the charge supercurrent in the system.
If you have Gnuplot installed, you can visualize the results by running `gnuplot` in the terminal, and then typing:

    :::gnuplot
    set xlabel 'Position'
    set ylabel 'Charge current'
    set yrange [-0.012:0.012]

    plot 'supercurrent.dat' using 1:2

The data files are formatted using tab-separated values (TSV files), where the first two columns correspond to the position and charge current, respectively.
If you wish to use a different visualization tool, it should be straight-forward to import the file `supercurrent.dat`, and plot the first two columns against each other.

Another quantity that might be interesting for this system, is to look at the superconducting phase as a function of position.
The file `correlation.dat` contains three data columns, which corresponds to the position, the gap magnitude, and the superconducting phase, respectively.
Thus, this information can be plotted as follows:

    :::gnuplot
    set xlabel 'Position'
    set ylabel 'Superconducting phase'
    set yrange [-0.25:0.25]

    plot 'correlation.dat' using 1:3


### Selfconsistent calculations
For the next example, we consider an S/F bilayer with spin-orbit coupling in the ferromagnet.
We take the superconductor length to be \( L_{\mathrm{S}} = 2.5\xi \), and treat its order parameter selfconsistently.
The ferromagnet has length \( L_{\mathrm{F}} = 0.5\xi \), magnetic exchange field \( \boldsymbol{m} = 3\Delta_0\boldsymbol{e}_x \), Rashba coupling \( \alpha\xi = 1 \), and Dresselhaus coupling \( \beta\xi = 1 \).
The spin-orbit coupling is in-plane; since the junction direction is along the \(z\)-axis, this corresponds to a Rashba Hamiltonian \( \mathcal{H}_{\mathrm{R}} \sim \alpha(p_y \sigma_x - p_x \sigma_y) \) and Dresselhaus Hamiltonian \( \mathcal{H}_{\mathrm{D}} \sim \beta(p_y\sigma_y - p_x\sigma_x) \).
The goal will be to calculate the local density of states in the junction.

An appropriate configuration file for this system would be:

    :::ini
    # S/F bilayer with SOC

    [superconductor]
      order:         2
      length:        2.5
      conductance_b: 0.3
    
    [ferromagnet]
      order:         1
      length:        0.5
      rashba:        1.0
      dresselhaus:   1.0
      magnetization: [3,0,0]
      conductance_a: 0.3

Save this file as e.g. `bilayer.conf`.
Note that there is no need to define boundary conditions for the "outer" interfaces: when no material is defined there, the program will automatically select vacuum boundary conditions.
The only thing needed to perform a selfconsistent calculation in the superconductor, is to specify a positive order for that layer (in this case 2).
The simulation can again be performed using the `converge` program:

    :::bash
    converge bilayer.conf &

The simulation program will alternate between solving the Usadel equation in the two layers.
Note that for the first few iterations, the order parameter in the superconductor is locked to \(\Delta = \Delta_0\).
This is because the state changes rapidly between iterations in the beginning, and selfconsistently updating the order parameter before the initial boundary conditions are approximately satisfied can slow down the convergence.
This "bootstrap procedure" continues until the change in the Riccati parameters between iterations is below 1%.
After that, the program start to perform proper selfconsistency iterations, where it both solves the Usadel equation and updates the order parameter.
After each such iteration, the physical observables in the system is written to output files, making it possible to preview the results.
Finally, every 8th iteration, a simple convergence acceleration procedure (Steffensen's method) is used predict what value the order parameter is converging towards, thus shortening the required simulation time.
Using IFort and a modern processor, the entire simulation should require roughly 30 min to converge (using a default error tolerance of \( 10^{-8} \) for the Riccati parameters).

Once the simulation has finished, the relevant output file is `density.dat`.
The first three columns of this file correspond to the position, energy, and the spin-independent local density of states, respectively.
The results are easily visualized as a contour plot in Gnuplot:

    :::gnuplot
    set xlabel  'Position'
    set ylabel  'Energy'
    set title   'Density of states'

    set yrange  [-3:3]
    set cbrange [0:2]

    set pm3d map
    splot 'density.dat' using 1:2:3



### Nonequilibrium calculations
@TODO This section is still under construction.

### Advanced functionality
@TODO This section is still under construction.

### Critical temperature
@TODO This section is still under construction.

### Phase diagrams
@TODO This section is still under construction.
