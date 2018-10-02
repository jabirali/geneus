title:  Tutorial
author: Jabir Ali Ouassou
date:   2018-09-03

### Introduction
In this tutorial, we demonstrate how to model some common superconducting spintronics devices using the GENEUS software suite.
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
In this section, we consider a voltage-biased N/S/N junction.
The voltage-biased contacts are used to inject a resistive charge current into the superconductor.
Once inside the superconductor, this resistive current decays as it is converted to a supercurrent.
Below, we intend to calculate and visualize that conversion process.
The superconductor is taken to have length \( L_{\mathrm{S}} = 5\xi \), and its order parameter has to be determined selfconsistently to obtain physically reasonable results.
Furthermore, due to the voltage bias, the system is necessarily out-of-equilibrium, which means that we also need to solve the kinetic equations numerically.
As for the normal metals, we will treat these as reservoirs at voltages \( eV/\Delta_0 = \pm 0.1 \).

An appropriate configuration file for the system above is:

    :::ini
    # Voltage-biased superconductor

    [conductor]
      order:           0
      voltage:        -0.1

    [superconductor]
      order:           1
      nonequilibrium:  T
      boost:           F
      length:          5.0
      conductance_a:   0.3
      conductance_b:   0.3

    [conductor]
      order:           0
      voltage:        +0.1

Save the configuration above as `voltage.conf`, and start the simulation:

    :::bash
    converge voltage.conf &

By default, a convergence acceleration method is invoked every 8th selfconsistency iteration.
For systems without phase gradients (i.e. charge supercurrents), this can speed up convergence by a factor 2-5x.
However, for systems with such phase gradients, the convergence acceleration procedure fails miserably, and should be disabled.
This is done by setting the parameter `boost` to false in the example above.

Note that you explicitly need to set `nonequilibrium` to true in each layer where you wish to solve the out-of-equilibrium kinetic equations.
Without this switch, the simulation program will assume that it is dealing with an equilibrium problem, and only solve the Usadel equation for the retarded propagators.

The example above is among the slowest single-layer junctions you can simulate with the code: using IFort and a modern processor, it takes roughly 90 min to converge completely.
This is because the code in general requires a large number of iterations to perform selfconsistent calculations with phase gradients, a problem that is further exacerbated by the lack of convergence acceleration.
If you wish to perform some faster numerical experiments to explore the out-of-equilibrium functionality, I would recommend testing N/N/N systems first.

Once the simulation finishes, the resistive current and supercurrent can be visualized together using Gnuplot:

    :::gnuplot 
    set xlabel 'Position'
    set ylabel 'Charge current'
    set yrange [0:0.002]

    plot 'supercurrent.dat' using 1:2 title 'Super', \
         'lossycurrent.dat' using 1:2 title 'Resistive'


### Advanced functionality
In this section, we introduce some more advanced capabilities of the configuration format: mathematical expressions and command-line arguments.
In the configuration file, you may replace any real number or vector by a mathematical expression.
The simulation program uses the [fparser](http://fparser.sourceforge.net/) library to parse these expressions, and supports all functionality of that library.
This basically means that you can use all elementary mathematical operators supported by Fortran itself (`**`, `*`, `/`, `+`, `-`), in addition to the most common elementary functions (`sin`, `cos`, `exp`, `log`, `sqrt`, etc.).
For more information, see the [fparser documentation](http://fparser.sourceforge.net/).

When defining mathematical expressions, there are two special variables available.
One of them is `pi`, which just refers to the mathematical constant (\( \pi = 3.1415\ldots \)).
The second is `z`, which refers to the position \( z \) inside a material.
The latter can be used to initialize physical fields using an analytical function of position, and is normalized so that that \( z = 0 \) at the left and \( z = 1 \) at the right interface of each material.

You can also feed the simulation programs command-line arguments.
These can then be referred to in the configuration file using the syntax `{1}` for the first command-line argument, `{2}` for the second, and so on.
In this way, the configuration files themselves specify how to interpret these arguments.
Note that if a command-line argument is referred to in the configuration file but not provided, the simulation program exits with an error.

Let us now consider a physical example where the functionality above might be useful.
Imagine an S/F/S Josephson junction, where the central layer is a helical ferromagnet like Ho.
Its magnetization texture can then be specified as a vector-valued function of position: \( \boldsymbol{m}(z) = [3\Delta_0 \cos(\pi z/2), 3\Delta_0 \sin(\pi z/2), \Delta_0] \).

As for the superconductors, these are modelled as reservoirs with fixed order parameters.
However, we wish to investigate the current-phase relation for multiple phase differences.
We therefore set the reservoir phases to \( \pm \delta\varphi/2 \), where the phase difference \( \delta\varphi \) is provided as a command-line argument `{1}`.
This allows us to perform simulations for different \( \delta\varphi \) in parallel without having to write multiple configuration files.

An appropriate configuration file is then:

    :::ini
    # Josephson junction with a helical ferromagnet.

    [superconductor]
      order:         0
      phase:        -{1}/2

    [ferromagnet]
      length:        3.0
      magnetization: [3*cos(pi*z/2), 3*sin(pi*z/2), 1]
      conductance_a: 0.3
      conductance_b: 0.3

    [superconductor]
      order:         0
      phase:        +{1}/2

Save this file as e.g. `helical.conf` and open a terminal in the same folder.
We can then use a Bash loop to start multiple parallel simulations with different \( \delta\varphi \in \{ 0.0, 0.1, \ldots, 1.0 \} \) in separate subdirectories:

    :::bash
    for n in $(seq 0.0 0.1 1.0); do
      mkdir sim_${n};
      cd sim_${n};
      converge ../helical.conf ${n} &
      cd ..;
    done

With IFort and an Intel Core i7, it takes roughly 3 min for all the simulations to complete.
Since the simulation results are scattered in multiple files in multiple folders, some postprocessing is required to collect the results afterwards.
The charge current is conserved, so we can extract it at any point in the ferromagnet; in this case, we use `tail -n 1` to extract it from the right end.
In Bash, these results can then be collected by running:

    :::bash
    for n in $(seq 0.0 0.1 1.0); do
      echo -n ${n} >> current.dat;
      tail -n 1 sim_${n}/supercurrent.dat >> current.dat
    done

Finally, the resulting current-phase relation can be plotted using:

    :::gnuplot
    set ylabel 'Charge current'
    set xlabel 'Phase difference'

    plot 'current.dat' using 1:3


### Critical temperature
In all of the examples above, we have focused on a single program `converge`, which essentially calculates the steady-state solution for a physical system with known parameters.
In this section, we will demonstrate how to calculate the critical temperature of a superconducting junction using a more specialized program `critical`.
The program uses a kind of binary search algorithm to drastically improve the calculation time needed to determine the critical temperature with a high precision.
The ideas behind this algorithm are explained in more detail in the appendix of [Scientific Reports 6, 29312 (2016)](https://www.nature.com/articles/srep29312).

As a model system, we consider an FI/S/FI spin-valve setup.
The ferromagnetic insulators are modelled as spin-active interfaces with spin-mixing conductances \( G_\varphi = 0.5G_{\mathrm{N}} \), where \( G_{\mathrm{N}} \) is the normal-state conductance of the superconductor.
Their magnetization directions are set to \( \boldsymbol{m} = [\cos(\pm\theta/2), \sin(\pm\theta/2), 0] \), where the relative magnetization angle \( \theta \) will be provided as a command-line argument.

An appropriate model for this spin-valve setup is:

    :::ini
    # FI/S/FI Spin-valve setup.

    [superconductor]
      # Material itself
      length:          1.0
    
      # Left interface
      magnetization_a: [cos(-{1}*pi/2),sin(-{1}*pi/2),0]
      spinmixing_a:    0.5
    
      # Right interface
      magnetization_b: [cos(+{1}*pi/2),sin(+{1}*pi/2),0]
      spinmixing_b:    0.5

Save this configuration file as e.g. `spinvalve.conf` and open a terminal in the same folder.
We then perform parallel simulations for magnetization angles \( \theta/\pi \in \{ 0.0, 0.1, \ldots, 1.0 \} \):

    :::bash
    for n in $(seq 0.0 0.1 1.0); do
      mkdir sim_${n};
      cd sim_${n};
      critical ../spinvalve.conf ${n} &
      cd ..;
    done

These simulations may take a few hours to complete.
When the simulations are finished, they will write both the command-line options used to invoke each process and the calculated critical temperature to `critical.dat` output files.
The critical temperature result \(T_{\mathrm{c}}\) is normalized to the critical temperature \(T_{\mathrm{cs}}\) of a bulk superconductor.
The results can be collected into a single output file as follows:

    :::bash
    cat sim_*/critical.dat | cut -f 2,3 > critical.dat

It is then trivial to plot the results in Gnuplot:

    :::gnuplot
    set xlabel 'Magnetization angle'
    set ylabel 'Critical temperature'

    plot 'critical.dat' u 1:2 notitle


### Phase diagrams
In this example, we demonstrate another specialized program `flow`, which can be used to rapidly map out the phase diagram of a superconducting system.
The program is based on analysing the "flow" of the order parameter during selfconsistency iterations; the concept is explained in more detail in the supplemental information of [arXiv:1803.07076](https://arxiv.org/abs/1803.07076).
The program will be illustrated using a spin-valve setup similar to the previous section.
However, in this case we will assume zero temperature and perpendicular magnetizations, and investigate how the spin-mixing conductance \( G_\varphi \) affects the stability of superconductivity.

An appropriate configuration file is:

    :::ini
    # FI/S/FI Spin-valve setup.

    [superconductor]
      # Material itself
      gap:             0.01
      length:          1.00
    
      # Left interface
      magnetization_a: [1,0,0]
      spinmixing_a:    {1}
    
      # Right interface
      magnetization_b: [0,1,0]
      spinmixing_b:    {1}

Save the file above as e.g. `spinvalve2.conf`, and run the following commands to start the simulations:

    :::bash
    for n in $(seq 0.0 0.1 1.0); do
      mkdir sim_${n};
      cd sim_${n};
      flow ../spinvalve2.conf ${n} &
      cd ..;
    done

According to the configuration file above, the order parameter is initially set to a small value \( \Delta = 0.01\Delta_0 \).
The simulation program then performs a fixed number of selfconsistency iterations, and determines whether the order parameter is spontaneously increasing (>1) or decreasing (<1) compared to its initial value.
This can be used to classify regions in parameter space as having a stable or unstable superconducting solution branch.

The simulations should take about 10 min to complete. 
The results can then be collected to a single output file:

    :::bash
    cat sim_*/flow.dat | cut -f 2,3 > flow.dat

They are then straight-forward to visualize in Gnuplot:

    :::gnuplot
    set xlabel 'Spin-mixing conductance'
    set ytics ('Normal' 0, 'Super' 1)

    plot 'flow.dat' u 1:($2 > 1) notitle

