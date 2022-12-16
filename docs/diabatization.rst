Diabatiziation
**************

Adiabatic states can be diabatized with a suitable rotation matrix, which in turn can be
determined by maximizing a cost funtion that depends on the rotation matrix and some selected
adiabatic properties (see below).
Pysisyphus implements property-based direct diabatizion, as outlined
by Subotnik [#subotnikBoys]_ and Truhlar et al. [#truhlarDQ]_

Currently, dipole moments, the trace of the primitive quadrupole tensor and
the electrostatic potential can be utilized for diabatization.

General Algorithm
-----------------

Diabatic states :math:`\{\Xi_A\}` can be obtained via unitary transformation
of adiabatic states :math:`\{\Psi_I\}`.

.. math::

   \Xi_A = \sum_{I=1} \Psi_I U_{IA}

which is equivalent to

.. math::

   \mathbf{\Xi} = \mathbf{U}^\intercal \Psi

with rotation matrix :math:`\mathbf{U}`, as well as adiabatic and diabatic
state vectors :math:`\mathbf{\Psi}` and :math:`\mathbf{\Xi}`.
Given two adiabatic states :math:`\Psi_1` and :math:`\Psi_2` :math:`\mathbf{U}`
is defined as

.. math::

   \mathbf{U} = \begin{bmatrix}
        \cos(\theta) & -\sin(\theta) \\
        \sin(\theta) & \cos(\theta) \\
    \end{bmatrix} ~ .

The related diabatic states :math:`\Xi_A` and :math:`\Xi_B` are

.. math::

   \begin{align}
     \mathbf{\Xi} &= \mathbf{U}^\intercal \mathbf{\Psi} \\
   \begin{bmatrix}
   \Xi_A \\
   \Xi_B \\
   \end{bmatrix} &= \begin{bmatrix}
        \cos(\theta) & \sin(\theta) \\
        -\sin(\theta) & \cos(\theta) \\
    \end{bmatrix}
   \begin{bmatrix}
   \Psi_1 \\
   \Psi_2 \\
   \end{bmatrix} \\
   \end{align}

.. math::
   \begin{align}
       \Xi_A &= \cos(\theta) \Psi_1 + \sin(\theta) \Psi_2 \\
       \Xi_B &= -\sin(\theta) \Psi_1 - \cos(\theta) \Psi_2 \\
   \end{align} ~ .

Similarly, diabatic expectation values are calculated as linear combination of
adiabatic expectation values:

.. math::
   \braket{\Xi_A | \hat{O} | \Xi_B} = \sum_{IJ} U_{IA} \braket{\Psi_I | \hat{O} | \Psi_J} U_{JB}

or in matrix notation

.. math::
   \mathbf{O}_\mathrm{dia} = \mathbf{U}^\intercal \mathbf{O}_\mathrm{adia} \mathbf{U} ~ .

By maximizing a cost function that depends on diabatic properties :math:`\mathbf{O}_\mathrm{dia}`,
a rotation matrix :math:`\mathbf{U}` suitable for diabatization can be determined.

With a known rotation matrix :math:`\mathbf{U}`, the diagonal adiabatic electronic energy
matrix :math:`\mathbf{V}` can be transformed to its diabatic equivalent :math:`\mathbf{D}`,
e.g, to determine diabatic couplings :math:`|D_{12}|`.

.. math::
   \begin{align}
     \mathbf{D} &= \mathbf{U}^\intercal \mathbf{V} \mathbf{U} \\
   \begin{bmatrix}
        D_{11} & D_{12} \\
        D_{12} & D_{22} \\
    \end{bmatrix}
    &= \begin{bmatrix}
        \cos(\theta) & \sin(\theta) \\
        -\sin(\theta) & \cos(\theta) \\
    \end{bmatrix}
   \begin{bmatrix}
        V_{11} & 0 \\
        0 & V_{22} \\
    \end{bmatrix}
   \begin{bmatrix}
        \cos(\theta) & -\sin(\theta) \\
        \sin(\theta) & \cos(\theta) \\
    \end{bmatrix}
   \end{align}

Depending on the employed properties, different cost functions :math:`f(\mathbf{U})` are used.
Subotnik [#subotnikBoys]_  proposed to maximize the magnitude of the diabatic dipole moments:

.. math::
   f_D(\mathbf{U}) = \sum_{A=1} |\braket{\Xi_A | \hat{\mu} | \Xi_A}|^2 ~ .

This is the same cost function as in classical Boys localization of molecular orbitals.
By also considering the trace of the primitive quadrupole tensor, :math:`f_D(\mathbf{U})`
is extended to :math:`f_{DQ}(\mathbf{U})`:

.. math::
   f_{DQ}(\mathbf{U}) = f_D(\mathbf{U})
    + \sum_{A=1} \sum_{j=1}^{N_Q} \alpha_j | \braket{\Xi_A | tr(\hat{Q}) | \Xi_A} |^2 ~ .

The contribution of the quadrupole moments to :math:`f(\mathbf{U})` is controlled by the
scaling factor :math:`\alpha_j`. Here, subscript :math:`j` refers to the j-th expansion
center. Currently, only one expansion center (:math:`N_Q = 1`) can be used in pysisyphus.
By default :math:`\alpha_j` is set to :math:`10 ~ a_0^{-2}`.

In some cases, considering only the dipole moments is not enough to discriminate different
adiabatic states and quadrupole moments have to be taken into account. Several examples are
outlined by Truhlar et al. [#truhlarDQ]_

Slight improvements may be possible by incorporating the electronic component of the
electrostatic potential (ESP) :math:`\Phi`.

.. math::
   f_{DQ\Phi}(\mathbf{U}) = f_{DQ}(\mathbf{U})
    + \sum_{A=1} \sum_{k=1}^{N_\Phi} \beta_k | \braket{\Xi_A | \hat{\Phi} | \Xi_A} |^2

Here again, :math:`N_\Phi` denotes the number of expansion centers where the ESP
is calculated. Similar to the quadrupole moments, only one expansion center is currently
supported in pysisyphus.
:math:`\beta_k` is a scaling factor that controls the contribution of the electrostatic
potential to the overall cost function :math:`f_{DQ\Phi}(\mathbf{U})`.

Maximization of cost function :math:`f` is achieved via repeated 2 x 2 rotations
(Jacobi sweeps) of diabatic state pairs. The formulas to calculate
the required rotation angle are given by Kleier [#kleierOpt]_ and Pulay et al. [#pulayOpt]_
The maximization is started with a unit-matrix guess for :math:`\mathbf{U}` and usually
converges quite rapidly.

A different approach based on unitary optimization is outlined
by Lehtola, [#lehtolaUnitary]_ but is not available in pysisyphus.


Diabatization is carried out via the `pysisdia` entry hook. See below for an example.

Example
-------

Diabatization of 3 adiabatic states using dipole moments requires the 3 adiabatic
permanent dipole moments (0,0; 1,1; 2,2), as well as the respective transition
dipole moments (0,1; 0,2; 1,2).

.. literalinclude :: ../tests/test_diabatization/guanine_indole_dia.yaml
   :language: yaml
   :caption:

Executing `pysisdia guanine_indole_dia.yaml` produces the following output:

.. code::

                              ###################
                              # D-DIABATIZATION #
                              ###################

    Dipole moments
    --------------
    [[[-1.253 -1.224  2.227]
      [-1.224 -1.691 -1.574]
      [ 2.227 -1.574  0.948]]

     [[ 0.27  -0.279  0.866]
      [-0.279 -0.653 -1.029]
      [ 0.866 -1.029  0.733]]

     [[-1.876 -0.325  0.735]
      [-0.325 -2.083 -0.517]
      [ 0.735 -0.517 -1.327]]]

    Starting Jacobi sweeps.
    000: P=169.63193103 dP=169.63193103
    001: P=172.49529007 dP=  2.86335904
    002: P=172.49529008 dP=  0.00000000
    Jacobi sweeps converged in 3 cycles.

    All energies are given in eV.

    Adiabatic energy matrix V
    -------------------------
    [[0.    0.    0.   ]
     [0.    0.717 0.   ]
     [0.    0.    0.905]]

    Rotation matrix U
    -----------------
    [[ 0.16289082  0.84920213  0.50231695]
     [-0.84058556  0.38601857 -0.38000734]
     [-0.51660672 -0.36034067  0.77670593]]
    det(U)=1.0000

    Diabatic energy matrix D = UᵀVU
    -------------------------------
    [[ 0.74814945 -0.06418359 -0.13410224]
     [-0.06418359  0.2243505  -0.35846691]
     [-0.13410224 -0.35846691  0.64950005]]

    Diabatic states Ξᵢ sorted by energy
    --------------------------------
    0: Ξ₁, 0.2244 eV
    1: Ξ₂, 0.6495 eV
    2: Ξ₀, 0.7481 eV

    Composition of diabatic states Ξᵢ
    ---------------------------------
    Ξ₀ = + 0.1629·Φ₀(GS) - 0.8406·Φ₁(LE) - 0.5166·Φ₂(CT)
    Ξ₁ = + 0.8492·Φ₀(GS) + 0.3860·Φ₁(LE) - 0.3603·Φ₂(CT)
    Ξ₂ = + 0.5023·Φ₀(GS) - 0.3800·Φ₁(LE) + 0.7767·Φ₂(CT)

    Weights U²
    ----------
    [[0.02653342 0.72114427 0.25232232]
     [0.70658408 0.14901034 0.14440558]
     [0.2668825  0.1298454  0.6032721 ]]

    Absolute diabatic couplings
    ---------------------------
    |D₀₁| = 0.0642 eV
    |D₀₂| = 0.1341 eV
    |D₁₂| = 0.3585 eV

YAML input
----------

All possible input options are outlined below.
The numbers are just dummy values.

.. code:: yaml

    # Adiabatic energies, of the states to diabatize.
    #
    # List of floats.
    # 	[V0, V1, ...]
    # May be absolute or relative energies.
    adiabatic_energies: [0.0, 0.6541, 0.7351]
    # Dipole moments.
    #
    # List of lists of length 5, each containing 2 integers followed by 3 floats.
    # 	[state0, state1, dpm_x, dpm_y, dpm_z]
    # The integers are adiabatic state indices, the 3 floats are the X, Y and Z
    # components of a dipole moment vector.
    # If both integers are the same, the 3 floats correspond to the permanent dipole moment
    # of a state. Otherwise, they correspond to a transition dipole moment
    # between two states.
    dipoles: [
     [0, 0, -1.0, -2.0, -3.0],
     [1, 1, -2.0, -1.0, -3.0],
     [0, 1, -2.0, -1.0, -3.0],
    ]
     
    #
    # Optional input below
    #

    # Adiabatic state labels.
    #
    # List of strings.
    # 	[label1, label2, ...]
    # Must be of same length as adiabatic energies.
    adiabiatc_labels: [GS, LE, CT]
    # Energy unit.
    #
    # Literal["eV", "Eh"]
    unit: eV
    # Trace of primitive quadrupole moment tensor.
    #
    # List of lists of length 3, each containing 2 integers followed by 1 float.
    # 	[state0, state1, tr_qpm] = [state0, state1, (qpm_xx + qpm_yy + qpm_zz) 
    # The same comments as for the dipole moments apply.
    quadrupoles: [
     [0, 0, 1.0],
     [1, 1, -14.0],
     [0, 1, 33.0],
    ]
    # Quadrupole moment scaling factor alpha in 1/a₀².
    #
    # Float
    alpha: 10.0
    # Electronic component of electrostatic potential in au
    epots: [
     [0, 0.48],
     [1, 0.12],
    ]
    # Electrostatic potential scaling factor beta in a₀.
    #
    # Float
    beta: 1.0


.. [#subotnikBoys] https://doi.org/10.1063/1.3042233
.. [#truhlarDQ] https://doi.org/10.1063/1.4894472, https://doi.org/10.1063/1.4948728
.. [#kleierOpt] https://doi.org/10.1063/1.1681683
.. [#pulayOpt] https://doi.org/10.1002/jcc.540140615
.. [#lehtolaUnitary] https://doi.org/10.1021/ct400793q
