# Readme

## 0.- Prerequisites

-- Install openmpi

```
sudo apt update
sudo apt upgrade
sudo apt install build-essential git gfortran openmpi-bin libopenmpi-dev
```

## 0.1- Installing

-- Download the repository from GitHub

```
git clone https://github.com/scsr-inpe/mpca-ann.git
```

-- Go to the folder and create the subdirectory build

```
cd mpcaann/
mkdir build
```

-- Make the project

```
make
```

## 1.- How to configure the experiments

The experiment configuration and the control parameters of the optimization algorithm are set in the file `configuration.ini` which is located in the subdirectory `config`:

### 1.1 Configuration of the neural networks model
* NCLASSES=8600, Number of patterns
* NCLASSESVALIDATION=1000, Number of patterns for validation
* NINPUTS=2, Number of inputs
* NOUTPUTS=1, Number of outputs
* TARGETERROR=1.0E-5, Value of error to be reached
* NEPOCHS=200, Maximum number of epochs (iterations of the training process)
* LOADWEIGHTSBIAS=1, Load predefined weights and biases.  If 0, weights and bias are fixed in 0.5. If 1, random weights and bias will be loaded. If 3, weights and bias are loaded from the file ``data/nn.initial``.
* HAVEVALIDATION=T, If T, train with validation; if F, do not use validation
* TRYINITIALARCHITECTURE=F, Try an initial (good?) architecture as initial solution. If F do not load. if V load configurations defined in 1.3.

### 1.2 Limits (search space) of the architecture of the neural network and the training parameters
* LOWER_HIDDEN_LAYERS=1,
* UPPER_HIDDEN_LAYERS=2,
* LOWER_FIRST_HIDDEN_LAYER=1,
* UPPER_FIRST_HIDDEN_LAYER=40,
* LOWER_SECOND_HIDDEN_LAYER=1,
* UPPER_SECOND_HIDDEN_LAYER=40,
* LOWER_ACTIVATION_FUNCTION=1,
* UPPER_ACTIVATION_FUNCTION=3,
* LOWER_ALPHA=1.0E-2,
* UPPER_ALPHA=0.9,
* LOWER_ETA=1.0E-2,
* UPPER_ETA=1.0,

### 1.3 Initial topology of the neural network
* INITIAL_HIDDEN_LAYERS=1, Number of hidden layers. Discrete [1; 2]
* INITIAL_FIRST_HIDDEN_LAYER=10, Number of neurons in the first hiden layer. Discrete  [1 - 40]
* INITIAL_SECOND_HIDDEN_LAYER=0, Number of neurons in the second hiden layer. Discrete  [1 - 40]
* INITIAL_ACTIVATION_FUNCTION=2, Activation function. Discrete [(1) logistic; (2) tangent; (3) gaussian]
* INITIAL_ALPHA=0.0, Momentum rate - alpha. Continuous [0.01 - 0.9]
* INITIAL_ETA=0.94, Learning rate - eta. Continuous [0.01 - 1]

### 1.4 Configuration of the optimization algorithm

Multi-particle collision algorithm

* VALUE_TO_REACH=1.0E-7, Stopping criterium. Minimum value of the objective function.
* PARTICLES_PROCESSOR=10, Number of particles in each processor
* MAXIMUM_NFE_MPCA=50000, Maximum number of function evaluation
* CYCLE_BLACKBOARD_MPCA=5000, Number of function evaluation to update the blackboard strategy. In blackboard, all particles know the best overall particle.
* NFE_EXPLOITATION_MPCA=1000, Number of function evaluation in the exploitation phase
* LOWER_EXPLOITATION_MPCA=0.7, Control parameter
* UPPER_EXPLOITATION_MPCA= 1.1, Control parameter.
* TYPE_PROBABILITY_MPCA=1, Probability distribution used in the scattering mechanism

Opposition (Hybrid variant)

* ENABLE_OPPOSITION=T, Enable opposition. If T, enables opposition. If F, dissables opposition
* TYPE_OPPOSITION="MPCA", Opposition type: MPCA (without opposition), ...

Hooke-Jeeve - Intensification phase (Hybrid variant)

* JUMPING_RATE_OPPOSITION=0.01, Probability to apply opposition in a particle
* EPSILON_HOOKE_JEEVES=1.0E-11, Control parameter
* RHO_HOOKE_JEEVES= 0.8, Control parameter
* MAXIMUM_NFE_HOOKE_JEEVES=1000, Maximum number of function evaluation

General

* VERBOSE=T, Displays extended information

## 2.- Data

### 2.1 Training data

The training data are loaded from the files ``data/x.txt`` and ``data/y.txt``, where ``x.txt`` represents the input set and ``y.txt`` the output set. Data is formated as F8.5 in Fortran, as ``      -0.07621``, with eight places before the decimal point and five decimal places.
Columns represent the inputs or the output, while rows contain the classes.

### 2.2 Cross-validation

The training dataset can be divided into two sets: training and cross-validation. The cross-validation dataset is used to check the overfitting problem and validate the trained network. Training should be stopped when the error in the validation starts to rise consistently.

The validation data are loaded from the files ``data/x_valid.txt`` and ``data/y_valid.txt``, where ``x_valid.txt`` represents the input set and ``y_valid.txt`` the output set. Data is formated as F8.5 in Fortran, as ``      -0.07621``, with eight places before the decimal point and five decimal places.
Columns represent the inputs or the output, while rows contain the classes.

## 3.- How to run the auto-configuration of the neural network

```
./runMPCA E P
```

where E is the total of experiments, and P is the number of processors to be used.

For example,

```
./runMPCA 5 4
```

runs the five experiments with four processor working.

## 4.- Output data

After the optimization process, the best configurations found will be available in the output folder.

The ``final.out`` shows a summary of the experiments. Each line represents an experiment:

```
 * Objective function minimum value;
 * Number of hidden layers;
 * Number of neurons in the first hiden layer;
 * Number of neurons in the second hiden layer;
 * Activation function. Discrete [(1) logistic; (2) tangent; (3) gaussian];
 * Momentum rate;
 * Learning rate.

  For example,

  8.540212E-02 1  5  0 1  3.771149E-01  3.906986E-01
```

The output folder contains files ``ann#.best``. This files contains the best configuration found in the experiment with number #.

An output file is formated as follows:

```
  8.536407E-02
  1
  7
wh1
  -10.86851    2.11830    0.97716    1.02200    1.62511   -0.56027    2.68275
  -10.36335    1.49072    1.27402    0.72460    1.27691   -0.41907    2.39177
bh1
    1.01121    1.19527    1.12059    1.19807    1.18954    0.52760   -0.76000
ws
  -10.56349
    2.14747
    1.07605
    0.80109
    1.53290
   -0.82369
    1.77848
bs
    4.72634
Alpha:        0.3566
Eta:        0.5252
Activation Function: LOGISTIC
MeanSquaredError:   2.788083E-03
MeanSquaredError - validation:   9.362167E-02
Epoch :       200
NFE:         22
```

## 5.- Using the best model found

### Generalization data

The generalization data are loaded from the files ``data/x_gen.txt`` and ``data/y_gen.txt``, where ``x_gen.txt`` represents the input set and ``y_gen.txt`` the output set. Data is formated as F8.5 in Fortran, as ``      -0.07621``, with eight places before the decimal point and five decimal places.
Columns represent the inputs or the output, while rows contain the classes.

The file ``config/annConfig.in`` has the configuration. the first line represent the classes, the second line is the number of inputs, and the third, the number of outputs

### Running the generalization process

To run the generalization simply execute the following command:

```
./annMLP E P
```

where E is the total of experiments and P is the number of processors to be used.

For example,

```
./annMLP 5 4
```

## Contributors

* Haroldo F. Campos Velho
* Juliana Anochi @juliana-anochi
* Sabrina Sambati @sabrinabms
* Eduardo Luz
* Reynier Hernández Torres @reynierhdez

## License

---
